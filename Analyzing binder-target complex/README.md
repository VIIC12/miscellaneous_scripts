# Analyzing Binder-Target Complex

This directory contains tools for analyzing protein-protein interactions, specifically focusing on binder-target complexes. The main script provides comprehensive structural analysis of protein complexes, calculating various metrics and properties.

## 🧬 Script Overview

The `analyse_binder_target.py` script performs detailed structural analysis of protein complexes, focusing on:
- Secondary structure analysis
- Surface area calculations
- Interface analysis
- Shape complementarity
- Glycan analysis
- Disulfide bond detection
- Clash detection
- And more...

## 📋 Features

### Binder Analysis
- Secondary structure composition (helix, sheet, loops)
- Sequence and length analysis
- Cysteine analysis (unpaired cysteines and disulfide bonds)
- Glycan analysis (surface vs. non-surface)
- Solvent accessible surface area (SASA) calculations
- Radius of gyration
- Disorder prediction

### Interface Analysis
- Interface residue identification
- Secondary structure composition at interface
- Shape complementarity scoring
- Contact analysis for specific residues

### Validation
- Clash detection between chains
- Interchain contact analysis for specified residues

## 🚀 Usage

```bash
python analyse_binder_target.py --pdb <input.pdb> --binder_chain <chain_id> --target_chains <chain_ids> [options]
```

### Required Arguments
- `--pdb`: Path to the input PDB file
- `--binder_chain`: Chain ID of the binder (e.g., 'A')
- `--target_chains`: Comma-separated list of target chains (e.g., 'H,L')

### Optional Arguments
- `--fixed_residues`: Space-separated list of residue numbers to check for contacts
- `--write_output`: Flag to write results to a JSON file

### Example
```bash
python analyse_binder_target.py --pdb complex.pdb --binder_chain A --target_chains H,L --fixed_residues "10 11 14 15" --write_output
```

## 📊 Output

The script outputs a comprehensive JSON structure containing:
- Binder properties (secondary structure, sequence, length, etc.)
- Interface properties (residues, secondary structure, shape complementarity)
- Validation results (contacts, clashes)

If `--write_output` is specified, results are saved to a JSON file named after the input PDB.

### Example Output Structure
```json
{
    "binder": {
        "helix": 0.0,
        "sheet": 0.0,
        "loops": 0.0,
        "sequence": "",
        "length": 0,
        "ss_seq": "",
        "condensed_structure": "",
        "intermediate_loops": 0,
        "unpaired_cys": 0,
        "disulfide": 0,
        "surface_glycans": 0,
        "non_surface_glycans": 0,
        "sasa_per_residue": 0.0,
        "fraction_disordered": 0.0,
        "radius_of_gyration": 0.0
    },
    "interface": {
        "helix": 0.0,
        "sheet": 0.0,
        "loops": 0.0,
        "binder_residues": [],
        "target_residues": [],
        "shape_complementarity": 0.0
    },
    "validation": {
        "interchain_contact": [],
        "clashes": false
    }
}
```

## 🔧 Dependencies

- Python 3.x
- PyRosetta
- BioPython
- NumPy
- FreeSASA

### Installation
```bash
conda create -n analyse_binder_target python=3.9 numpy freesasa && conda activate analyse_binder_target && pip install pyrosetta && pip install Bio
```

## 📝 Notes

- The script requires PyRosetta to be properly installed and initialized
- For glycan analysis, the script looks for NXT and NXS motifs
- Surface accessibility is calculated using FreeSASA
- Shape complementarity is calculated using PyRosetta's ShapeComplementarityFilter


## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📧 Contact

Tom U. Schlegel - Institute for Drug Discovery, Leipzig University 