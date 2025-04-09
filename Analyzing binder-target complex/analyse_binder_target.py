#!/usr/bin/env python3
"""
Author:     Tom U. Schlegel
Date:       2025-04-09
Name:       structure_analyse
Info:       Analyse structure and get various scores.
"""

import os
import sys
import argparse
import json

import numpy as np
import freesasa

import pyrosetta
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.core.scoring.sasa import SasaCalc
from pyrosetta.rosetta.core.pose import get_center_of_mass, radius_of_gyration, Pose
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, InterGroupInterfaceByVectorSelector, TrueResidueSelector, OrResidueSelector
)
from pyrosetta.rosetta.protocols.moves import DsspMover
from pyrosetta.rosetta.protocols.simple_filters import (ShapeComplementarityFilter)

from Bio.PDB import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

pyrosetta.init("-mute all", silent=True)

# Reference values from Wilke  et al. (2013) - Empirical
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5314839/#:~:text=The%20relative%20solvent%20accessible%20surface%20area%20(rASA)%20of%20a%20residue,a%20solvent%20molecule%20%5B22%5D
MAX_ASA_REF = {
    'ALA': 121.0,
    'ARG': 265.0,
    'ASN': 187.0,
    'ASP': 187.0,
    'CYS': 148.0,
    'GLU': 214.0,
    'GLN': 214.0,
    'GLY': 97.0,
    'HIS': 216.0,
    'ILE': 195.0,
    'LEU': 191.0,
    'LYS': 230.0,
    'MET': 203.0,
    'PHE': 228.0,
    'PRO': 154.0,
    'SER': 143.0,
    'THR': 163.0,
    'TRP': 264.0,
    'TYR': 255.0,
    'VAL': 165.0,
}

def secondary_structure_binder(pose: pyrosetta.Pose, chain_id: str) -> tuple[float, float, float, str, int, str, str, int]:
    """
    Calculate the secondary structure composition of a chain from a pdbfile.
    
    Args:
        pose: PyRosetta pose object
        chain_id: Chain identifier
        
    Returns:
        SecondaryStructureAnalysis with the following fields:
        - frac_helix: Fraction of helical structure
        - frac_sheet: Fraction of beta sheet structure
        - frac_loops: Fraction of loop structure
        - binder_seq: Full sequence of the chain
        - binder_length: Total length of the chain
        - secondary_structure: Full secondary structure sequence
        - condensed_structure: Simplified secondary structure sequence
        - intermediate_loops_count: Number of loops between secondary structure elements
    """

    try:
        Dssp = pyrosetta.rosetta.protocols.moves.DsspMover()
        Dssp.apply(pose)
        
        # Get sequence of specified chain
        binder_seq = ""
        for i in range(1, pose.size() + 1):
            if pose.pdb_info().chain(i) == chain_id:
                binder_seq += pose.residue(i).name1()

        # Get secondary structure in one pass
        ss = pose.secstruct()
        ss_seq = []
        for i in range(1, pose.size() + 1):
            if pose.pdb_info().chain(i) == chain_id:
                ss_seq.append(ss[i-1])
        ss_seq = ''.join(ss_seq)

        # Calculate fractions in one pass
        total = len(ss_seq)
        counts = {'H': 0, 'E': 0, 'L': 0}
        condensed_structure = []
        prev_ss = None
        
        for ss_type in ss_seq:
            counts[ss_type] += 1
            if ss_type != prev_ss:
                condensed_structure.append(ss_type)
                prev_ss = ss_type

        # Count intermediate loops in one pass
        loop_count = sum(
            1 for i in range(1, len(condensed_structure) - 1)
            if (condensed_structure[i] == 'L' and 
                condensed_structure[i-1] in 'EH' and 
                condensed_structure[i+1] in 'EH')
        )

        frac_helix = counts['H'] / total
        frac_sheet = counts['E'] / total
        frac_loops = counts['L'] / total
        
        binder_length = total
        ss_seq = ''.join(ss_seq)
        condensed_structure = ''.join(condensed_structure)
        intermediate_loops_count = loop_count

        return frac_helix, frac_sheet, frac_loops, binder_seq, binder_length, ss_seq, condensed_structure, intermediate_loops_count

    
    except Exception as e:
        print(f"Error in calc_secondary_structure: {str(e)}")
        return None

def calculate_distance(residue1, residue2) -> float:
    """
    Calculate distance between SG atoms of two cysteine residues

    Args:
        residue1: BioPython residue
        residue2: BioPython residue

    Returns:
        Distance between SG atoms of two cysteine residues
    """
    try:
        return residue1['SG'] - residue2['SG']
    except KeyError:
        return float('inf')

def analyze_cysteines(structure: PDBParser, chain_id: str) -> tuple[int, int]:
    """
    Analyze cysteines in the structure to identify disulfide bonds and unpaired cysteines.
    Returns tuple of (unpaired_cysteines_count, disulfide_bonds_count)

    Args:
        structure: BioPython structure
        chain_id: Chain ID

    Returns:
        Tuple of (unpaired_cysteines_count, disulfide_bonds_count)
    """
    # Distance threshold for disulfide bonds (in Angstroms)
    DISULFIDE_THRESHOLD = 2.2

    # Get all cysteine residues in the specified chain
    model = structure[0]
    chain = model[chain_id]
    cysteines = [res for res in chain if res.get_resname() == "CYS"]
    
    if not cysteines:
        return 0, 0  # No cysteines found
    
    # Track which cysteines are involved in disulfide bonds
    involved_in_bonds = set()
    disulfide_count = 0
    
    # Check all pairs of cysteines
    for i, cys1 in enumerate(cysteines):
        for cys2 in cysteines[i+1:]:
            distance = calculate_distance(cys1, cys2)
            if distance <= DISULFIDE_THRESHOLD:
                involved_in_bonds.add(cys1.get_id()[1])  # Get residue number
                involved_in_bonds.add(cys2.get_id()[1])
                disulfide_count += 1
    
    # Count unpaired cysteines
    unpaired_count = len(cysteines) - len(involved_in_bonds)
    
    return unpaired_count, disulfide_count

def calc_glycans(pose: pyrosetta.Pose, chain_id: str) -> tuple[int, int]:
    """
    Calculate the number of surface-accessible and non-surface-accessible glycans in a protein structure.

    Args:
        pose: PyRosetta pose object
        chain_id: Chain ID

    Returns:
        Tuple of (surface_glycans, non_surface_glycans)
    """

    sasa_calculator = SasaCalc()
    sasa_calculator.calculate(pose)
    sasa_values = sasa_calculator.get_residue_sasa()

    sasa_cutoff = 10.0 #! Adjust as needed
    surface_residues = [i for i, sasa in enumerate(sasa_values, start=1) if sasa > sasa_cutoff]

    # Get residues for specified chain or all chains
    if chain_id:
        chain_num = ord(chain_id) - ord('A') + 1  # Convert chain letter to number (A=1, B=2, etc)
        chain_begin = pose.conformation().chain_begin(chain_num)
        chain_end = pose.conformation().chain_end(chain_num)
        residue_range = range(chain_begin, chain_end + 1)
    else:
        residue_range = range(1, pose.size() + 1)

    asn_residues = [i for i in residue_range if pose.residue(i).name1() == "N"]

    nxt_surface_motifs = []
    nxs_surface_motifs = []
    nxt_non_surface_motifs = []
    nxs_non_surface_motifs = []

    for res in asn_residues:
        if res < pose.size() - 1:
            if pose.residue(res + 2).name1() == "T":
                if res in surface_residues:
                    nxt_surface_motifs.append(res)
                else:
                    nxt_non_surface_motifs.append(res)

            if pose.residue(res + 2).name1() == "S":
                if res in surface_residues:
                    nxs_surface_motifs.append(res)
                else:
                    nxs_non_surface_motifs.append(res)

    surface_glycans=len(nxt_surface_motifs) + len(nxs_surface_motifs)
    non_surface_glycans=len(nxt_non_surface_motifs) + len(nxs_non_surface_motifs)
    
    return surface_glycans, non_surface_glycans

def calculate_surface_metrics(pdb_path: str, chain_id) -> tuple[float, np.ndarray]:
    """
    Calculate the surface area of a protein structure.

    Args:
        pdb_path: Path to the PDB file
        chain_id: Chain ID

    Returns:
        Tuple of (sasa_per_residue, rasa_array)
    """
    #Set up the classifier and parameters
    classifier = freesasa.Classifier.getStandardClassifier('naccess')
    parameters = freesasa.Parameters({'algorithm': freesasa.LeeRichards, 'n-slices': 100})
    freesasa.setVerbosity(1)
    
    # Calculate the SASA, return result object
    sasa_calc_results = freesasa.calc(freesasa.Structure(pdb_path, classifier), parameters)

    # Get total SASA
    total_sasa = sasa_calc_results.totalArea()
    #print(f"total_sasa: {total_sasa}")

    # Get SASA of specified chain only
    sasa_subselection = freesasa.selectArea([f'chain, Chain {chain_id}'], freesasa.Structure(pdb_path), sasa_calc_results)
    #print(f"sasa_subselection['chain']: {sasa_subselection['chain']}")

    # Get SASA per residue of chain
    num_residues = len(sasa_calc_results.residueAreas()[chain_id].keys())
    #print(f"num_residues: {num_residues}")

    sasa_per_residue = sasa_subselection['chain'] / num_residues
    #print(f"sasa_per_residue: {sasa_per_residue}")

    # RASA
    rasa_values = []
    residue_numbers = list(sasa_calc_results.residueAreas()[chain_id].keys())
    for residue_num in residue_numbers:
        asa = sasa_calc_results.residueAreas()[chain_id][residue_num]
        #print(f"Residue {residue_num}: {asa.residueType} {asa.total}")
        
        max_asa = MAX_ASA_REF[asa.residueType]

        rasa = asa.total / max_asa if max_asa > 0 else 0
        #print(f"RASA: {rasa}")
        rasa_values.append(rasa)

    rasa_array = np.array(rasa_values)

    return sasa_per_residue, rasa_array

def estimate_fraction_disordered(rasa_array, threshold=0.2) -> float:
    """
    Estimate the fraction of disordered residues in a protein structure.
    threshold of 0.2 (or 20% relative accessibility) is an empirically determined value.
    based on studies that have found this level of exposure to be a good indicator of potential disorder.
    
    Args:
        rasa_array: Array of relative accessible surface areas
        threshold: Threshold for disordered residues

    Returns:
        Fraction of disordered residues
    """
    disordered = rasa_array > threshold
    return float(np.mean(disordered))

def get_interface_residues(pose: pyrosetta.Pose, chain_binder: str, chain_targets: list[str], cb_dist: float=7.0) -> list[str]:
    """
    Identify interface residues between binder and target(s).
    
    Args:
        pose: PyRosetta pose object
        chain_binder: Chain ID of the binder
        chain_targets: List of target chain IDs
        cb_dist: Distance cutoff for interface residues in Angstroms

    Returns:
        Tuple of two lists:
        - List of target residues in the format "chainIDresidueNumber"
        - List of binder residues in the format "chainIDresidueNumber"
    """

    interface = InterGroupInterfaceByVectorSelector()
    
    # Create selector for all target chains
    target_selector = ChainSelector(chain_targets[0])
    for chain in chain_targets[1:]:
        additional_chain = ChainSelector(chain)
        combined_selector = OrResidueSelector(target_selector, additional_chain)
        target_selector = combined_selector
        
    interface.group1_selector(target_selector)
    interface.group2_selector(ChainSelector(chain_binder))
    interface.cb_dist_cut(cb_dist)
    vec = interface.apply(pose)
    
    binder_residues = []
    target_residues = []
    for i, pos in enumerate(vec):
        if pos:
            try:
                pdb_info = pose.pdb_info().pose2pdb(i + 1).split()
                residue_num = pdb_info[0]
                chain = pdb_info[1]
                residue_str = f"{chain}{residue_num}"
                if chain == chain_binder:
                    binder_residues.append(residue_str)
                elif chain in chain_targets:
                    target_residues.append(residue_str)
            except Exception as e:
                print(f"Warning: Could not process residue at position {i + 1}: {e}")

    return target_residues, binder_residues

def get_interface_secondary_structure_fractions(pose: pyrosetta.Pose, chain_binder: str, chain_targets: list[str], cb_dist: float=7.0) -> dict:
    """
    Calculate secondary structure fractions at the interface for chain_binder.

    Args:
        pose: PyRosetta pose object
        chain_binder: Chain ID of the binder
        chain_targets: List of target chain IDs
        cb_dist: Distance cutoff for interface residues in Angstroms

    Returns:
        Dictionary containing the secondary structure fractions
    """
    interface_residues, binder_residues = get_interface_residues(pose, chain_binder, chain_targets, cb_dist)
    interface_residues_chain_binder = [res for res in interface_residues if res.startswith(chain_binder)]    
    DsspMover().apply(pose)
    ss = pose.secstruct()
    interface_ss = []

    for res in interface_residues_chain_binder:
        try:
            res_num = int(''.join(filter(str.isdigit, res)))
            pose_idx = pose.pdb_info().pdb2pose(chain_binder, res_num)
            if pose_idx > 0:  # Valid pose index
                interface_ss.append(ss[pose_idx - 1])
        except Exception as e:
            print(f"Warning: Could not process residue '{res}': {e}")
            continue

    total = len(interface_ss)
    if total == 0:
        return 0.0, 0.0, 0.0
    
    frac_helix_interface = interface_ss.count('H') / total
    frac_sheet_interface = interface_ss.count('E') / total
    frac_loops_interface = interface_ss.count('L') / total

    return frac_helix_interface, frac_sheet_interface, frac_loops_interface

def calc_shape_complementarity(pose: pyrosetta.Pose, chain_binder: str, chain_targets: list[str]) -> float:
    """
    Calculate shape complementarity between chains.

    Args:
        pose: PyRosetta pose object
        chain_binder: Chain ID of the binder
        chain_targets: List of target chain IDs

    Returns:
        Shape complementarity score
    """
    shape = ShapeComplementarityFilter()
    target_selector = ChainSelector(chain_targets[0])
    for chain in chain_targets[1:]:
        additional_chain = ChainSelector(chain)
        combined_selector = OrResidueSelector(target_selector, additional_chain)
        target_selector = combined_selector
        
    shape.selector1(target_selector)
    shape.selector2(ChainSelector(chain_binder))
    return shape.score(pose)

def get_radius_of_gyration(pose: pyrosetta.Pose, chain: str) -> float:
    """
    Calculate the radius of gyration for the specified chain.

    Args:
        pose: PyRosetta pose object
        chain: Chain ID

    Returns:
        Radius of gyration
    """
    chain_selector = ChainSelector(chain)
    selected_residues = chain_selector.apply(pose)
    chain_pose = Pose()
    for i in range(1, pose.total_residue() + 1):
        if selected_residues[i]:
            chain_pose.append_residue_by_bond(pose.residue(i))
    selector = TrueResidueSelector()
    selected_residues_new_pose = selector.apply(chain_pose)
    center = get_center_of_mass(chain_pose)
    rog = radius_of_gyration(chain_pose, center, selected_residues_new_pose)

    return rog

def check_clashes(pose_biopython, binder_chain: str, target_chains: list[str], cutoff: float=1.5) -> bool:
    """
    Checks for steric clashes between the binder and target chains using BioPython.
    
    Args:
        pose_biopython: BioPython structure
        binder_chain: Chain ID of the binder
        target_chains: List of target chain IDs
        cutoff: Distance cutoff for clashes in Angstroms

    Returns:
        True if there is a clash (atoms closer than cutoff), False otherwise
    """
    try:        
        # Get binder atoms
        binder_atoms = list(pose_biopython[0][binder_chain].get_atoms())
        
        # Get all target chain atoms
        target_atoms = []
        for chain in target_chains:
            target_atoms.extend(list(pose_biopython[0][chain].get_atoms()))
        
        # Use neighbor search to find clashes
        ns = NeighborSearch(target_atoms)
        for atom in binder_atoms:
            # Find all atoms within cutoff distance
            neighbors = ns.search(atom.get_coord(), cutoff)
            if neighbors:
                return True  # Clash found
                
        return False
        
    except Exception as e:
        print(f"Error in remove_clashes: {str(e)}")
        return False

def check_interchain_contact(pose_biopython, binder_chain: str, target_chains: list[str], fixed_residues_string: str, cutoff: float=3) -> list[bool]:
    """
    Checks if specific residues in the binder chain have contacts with any target chain atoms.
    
    Args:
        pose_biopython: BioPython structure
        binder_chain: Chain ID of the binder
        target_chains: List of target chain IDs
        fixed_residue: Space-separated string of residue numbers (e.g. "10 11 14 15")
        cutoff: Distance cutoff for contacts in Angstroms

    Returns:
        List of booleans indicating contact status for each residue in order
        e.g., [True, False, True] for "11 14 15"
    """
    try:
        # Convert fixed_residue string to list of integers
        fixed_residues = [int(res) for res in fixed_residues_string.split()]
        
        # Initialize results list
        contact_results = []
            
        # Get all target chain atoms
        target_atoms = []
        for chain in target_chains:
            target_atoms.extend(list(pose_biopython[0][chain].get_atoms()))
        
        ns = NeighborSearch(target_atoms)
        
        # Check each fixed residue in binder chain
        for res_num in fixed_residues:
            try:
                binder_residue = None
                for residue in pose_biopython[0][binder_chain]:
                    if residue.id[1] == res_num:
                        binder_residue = residue
                        break
                        
                if not binder_residue or 'CA' not in binder_residue:
                    contact_results.append(False)
                    continue
                    
                binder_ca = binder_residue['CA']
                
                # Use neighbor search to find contacts
                contacts = ns.search(binder_ca.get_coord(), cutoff)
                contact_results.append(len(contacts) > 0)
                    
            except Exception as e:
                print(f"Warning: Could not process residue {res_num}: {e}")
                contact_results.append(False)
                
        return contact_results
        
    except Exception as e:
        print(f"Error in check_interchain_contact: {str(e)}")
        return [False] * len(fixed_residues_string.split())

def main():
    parser = argparse.ArgumentParser(description="Filter PDBs based on various properties.")
    parser.add_argument('--pdb', required=True, help="Path to a specific PDB file to process", type=str)
    parser.add_argument('--binder_chain', required=True, help="Chain ID of the binder (e.g. A)", type=str)
    parser.add_argument('--target_chains', required=True, help="Comma-separated list of target chains (e.g. H,L for antibody)", type=str)
    parser.add_argument('--fixed_residues', required=False, help="Residue number of the fixed residue (e.g. '10 11 14 15\')", type=str)
    parser.add_argument('--delimiter', required=False, help="Delimiter for the fixed residues (default: ,)", type=str, default=",")
    parser.add_argument('--write_output', required=False, action='store_true', help="Write output to PDB-Name.json file")
    parser.epilog = "Example: python structure_analyse.py --pdb input.pdb --binder_chain A --target_chains H,L --fixed_residues '10 11 14 15'"
    args = parser.parse_args()

    # Split target chains into a list
    target_chains = args.target_chains.split(',')

    # BioPython
    parser = PDBParser(QUIET=True)
    pose_biopython = parser.get_structure('protein', args.pdb)
    
        # Validate chains exist in structure
    available_chains = {chain.id for chain in pose_biopython[0]}
    
    if args.binder_chain not in available_chains:
        print(f"Error: Binder chain '{args.binder_chain}' not found in PDB. Available chains: {', '.join(sorted(available_chains))}")
        sys.exit(1)
        
    missing_targets = [chain for chain in target_chains if chain not in available_chains]
    if missing_targets:
        print(f"Error: Target chain(s) {', '.join(missing_targets)} not found in PDB. Available chains: {', '.join(sorted(available_chains))}")
        sys.exit(1)

    # PyRosetta
    pose_pyrosetta = pyrosetta.pose_from_pdb(args.pdb)

    # Analyse binder
    frac_helix_binder, frac_sheet_binder, frac_loops_binder, binder_seq, binder_length, ss_seq_binder, condensed_structure_binder, intermediate_loops_count_binder = secondary_structure_binder(pose_pyrosetta, args.binder_chain)
    unpaired_cys, disulfide = analyze_cysteines(pose_biopython, args.binder_chain)
    surface_glycans, non_surface_glycans = calc_glycans(pose_pyrosetta, args.binder_chain)
    sasa_per_residue, rasa_array = calculate_surface_metrics(args.pdb, args.binder_chain)
    fraction_disordered = estimate_fraction_disordered(rasa_array)
    radius_of_gyration = get_radius_of_gyration(pose_pyrosetta, args.binder_chain)

    # Analyse interface
    target_residues, binder_residues = get_interface_residues(pose_pyrosetta, args.binder_chain, target_chains)
    frac_helix_interface, frac_sheet_interface, frac_loops_interface = get_interface_secondary_structure_fractions(pose_pyrosetta, args.binder_chain, target_chains)
    shape_complementarity = calc_shape_complementarity(pose_pyrosetta, args.binder_chain, target_chains)

    # Control
    interchain_contact = None
    if args.fixed_residues:
        interchain_contact = check_interchain_contact(pose_biopython, args.binder_chain, target_chains, args.fixed_residues)
    clashes = check_clashes(pose_biopython, args.binder_chain, target_chains)

    # Results
    results = {
        "binder": {
            "helix": round(frac_helix_binder, 4),
            "sheet": round(frac_sheet_binder, 4),
            "loops": round(frac_loops_binder, 4),
            "sequence": binder_seq,
            "length": binder_length,
            "ss_seq": ss_seq_binder,
            "condensed_structure": condensed_structure_binder,
            "intermediate_loops": intermediate_loops_count_binder,
            "unpaired_cys": unpaired_cys,
            "disulfide": disulfide,
            "surface_glycans": surface_glycans,
            "non_surface_glycans": non_surface_glycans,
            "sasa_per_residue": round(sasa_per_residue, 4),
            "fraction_disordered": round(fraction_disordered, 4),
            "radius_of_gyration": round(radius_of_gyration, 4),
        },
        "interface": {
            "helix": round(frac_helix_interface, 4),
            "sheet": round(frac_sheet_interface, 4),
            "loops": round(frac_loops_interface, 4),
            "binder_residues": f"{args.delimiter}".join(binder_residues),
            "target_residues": f"{args.delimiter}".join(target_residues),
            "shape_complementarity": round(shape_complementarity, 4)
        },
        "validation": {
            "interchain_contact": f"{args.delimiter}".join(str(x) for x in interchain_contact),
            "clashes": clashes
        }
    }

    if args.write_output:
        # Get PDB filename without path and extension
        pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
        output_file = f"{pdb_name}.json"
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=4)
    else:
        print(json.dumps(results))
    
if __name__ == "__main__":
    main()