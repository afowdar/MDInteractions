# MDInteractions: A python package for the distance-based analysis of intra- and inter-protein interactions.

MDInteractions is a Python package designed to analyze atomic-level intra- and inter-protein interactions from molecular dynamics (MD) simulations. It allows users to specify interaction parameters such as distance thresholds to detect and track contacts between atoms or residues within and between proteins over time.
 

## Installation

To install the package, run:

```
pip install MDInteractions
```

## Modules

- **protein_interactions**: The protein_interactions module investigates interactions within a protein group or between two protein groups, identifying residues that are within a user-defined atom-distance threshold at each frame. The final pairs recorded in the CSV file are those that consistently remained within the specified cutoff distance across all of the processed frames.
- **mean_distance**: The mean_distance module calculates the average distance between all possible pairs of user-defined atoms across one group for intra-protein interactions and across two groups for inter-protein interactions, throughout the user-defined frames of the trajectory. 

## Example Usage

```python
# Intra-protein interactions: Analyze contacts within a single group
from MDInteractions import protein_interactions

intra_interactions = protein_interactions(
    gro_file="md.gro",                      # Input structure file (.gro)
    xtc_file="md.xtc",                      # Input trajectory file (.xtc)
    ndx_file="index.ndx",                   # GROMACS index file defining atom groups
    start_frame=8000,                       # Frame to start analysis from
    end_frame=9000,                         # Frame to end analysis at
    group_ID=10,                            # Single group ID (intra-protein mode)
    group1_atom_name="CB",                  # One or more atom names as comma-separated string (e.g., 'CB' or 'CB, CA') for first atom set 
    group2_atom_name="CB",                  # Same as above, for second atom set 
    cutoff=7,                             # Distance cutoff for defining an interaction
    residue_specific_atoms=[("GLY", "CA")], # list of residue-specific atom tuples, e.g. '[ ("GLY", "CA"), ("ALA", "CB") ]'
    give_res_name=True,                     # Include residue names in the output
    give_atom_name=True                     # Include atom names in the output
)

intra_interactions.analyze()                # Run the interaction analysis
```

```python
# Inter-protein interactions: Analyze contacts between two distinct groups
from MDInteractions import protein_interactions

inter_interactions = protein_interactions(
    gro_file="md.gro",                    
    xtc_file="md.xtc",                    
    ndx_file="index.ndx",                 
    start_frame=8000,                     
    end_frame=9000,                       
    group1_ID=10,                           # First group ID for inter-protein 
    group2_ID=11,                           # Second group ID for inter-protein 
    group1_atom_name="CB",                  # One or more atom names as comma-separated string (e.g., 'CB' or 'CB, CA')
    group2_atom_name="CB",                  # Same as above, for group 2
    cutoff=15,                              # Larger cutoff for inter-group interactions
    residue_specific_atoms=[("GLY", "CA")], # list of residue-specific atom tuples, e.g. '[ ("GLY", "CA"), ("ALA", "CB") ]'
    give_res_name=True,                 
    give_atom_name=True
)

inter_interactions.analyze()                # Run the interaction analysis
```

```python
# Intra-protein mean distance: Average distance of selected atoms within a group 
from MDInteractions import mean_distance

intra_meandistance = mean_distance(
    gro_file="md.gro",                    
    xtc_file="md.xtc",                    
    ndx_file="index.ndx",                 
    start_frame=8000,                     
    end_frame=9000,                       
    group_ID=10,                            # Single group ID (intra-protein)
    group1_atom_name="CB",                  # One or more atom names as comma-separated string (e.g., 'CB' or 'CB, CA') for first atom set within group
    group2_atom_name="CB",                  # Same as above, for second atom set within group
    residue_specific_atoms=[("GLY", "CA")], # Python-style list of tuples. E.g., '[("GLY", "CA"), ("LYS", "NZ")]'
    give_res_name=True,                     # Include residue names in output
    give_atom_name=True,                    # Include atom names in output

)

intra_meandistance.analyze()                # Compute average distances
```

```python
# Inter-protein mean distance: Average distance between selected atoms across two groups
from MDInteractions import mean_istance

inter_meandistance = mean_distance(
    gro_file="md.gro",                    
    xtc_file="md.xtc",                    
    ndx_file="index.ndx",                 
    start_frame=8000,                     
    end_frame=9000,                       
    group1_ID=10,                           # First group ID for inter-protein mode
    group2_ID=11,                           # Second group ID for inter-protein mode
    group1_atom_name="CB",                  # One or more atom names as comma-separated string (e.g., 'CB' or 'CB, CA')
    group2_atom_name="CB",                  # Same as above, for group 2
    residue_specific_atoms=[("GLY", "CA")], # Python-style list of tuples. E.g., '[("GLY", "CA"), ("LYS", "NZ")]'
    give_res_name=True,                   
    give_atom_name=True,                  
    )

inter_meandistance.analyze()               # Compute average distances
```

## Command Line Interface 
```bash
# Intra-protein interaction analysis (within one group)

# Input GROMACS .gro structure file
protein_interactions \
  --gro_file "md.gro" \
# Input GROMACS .xtc trajectory file
  --xtc_file "md.xtc" \
# GROMACS index file defining groups
  --ndx_file "index.ndx" \
# Starting frame index for analysis
  --start_frame 8000 \
# Ending frame index for analysis
  --end_frame 9000 \
# Single group ID (intra-protein)
  --group_ID 10 \
# Atom name in group (e.g., beta carbon)
  --group1_atom_name CB \
# Atom name for distance comparison (same group)
  --group2_atom_name CB \
# Distance cutoff in Ångströms
  --cutoff 7.5 \
# JSON list of residue-specific atom tuples, e.g. '[ ["GLY", "CA", "CB"], ["ALA", "CB"] ]'
  --residue_specific_atoms '[["GLY", "CA"]]' \
# Include residue names in output
  --give_res_name=True \
# Include atom names in output
  --give_atom_name=True
```

```bash
# Inter-protein interaction analysis

# Input GROMACS .gro structure file
protein_interactions \
  --gro_file "md.gro" \
# Input GROMACS .xtc trajectory file
  --xtc_file "md.xtc" \
# GROMACS index file defining groups
  --ndx_file "index.ndx" \
# Starting frame index for analysis
  --start_frame 8000 \
# Ending frame index for analysis
  --end_frame 9000 \
# First group ID for inter-protein
  --group1_ID 10 \
# Second group ID for inter-protein
  --group2_ID 11 \
# Atom name for group 1
  --group1_atom_name CB \
# Atom name for group 2
  --group2_atom_name CB \
# Larger cutoff for inter-group distances
  --cutoff 15 \
# JSON list of residue-specific atom tuples, e.g. '[ ["GLY", "CA"], ["ALA", "CB"] ]'
  --residue_specific_atoms '[["GLY", "CA"]]' \
# Include residue names in output
  --give_res_name=True \
# Include atom names in output
  --give_atom_name=True
```

```bash
# Mean distance calculation for intra-protein 

# Input GROMACS .gro structure file
mean_distance \
  --gro_file "md.gro" \
# Input GROMACS .xtc trajectory file
  --xtc_file "md.xtc" \
# GROMACS index file defining groups
  --ndx_file "index.ndx" \
# Starting frame index for analysis
  --start_frame 8000 \
# Ending frame index for analysis
  --end_frame 9000 \
# Single group ID (intra-protein)
  --group_ID 10 \
# Atom name in group (e.g., beta carbon)
  --group1_atom_name CB \
# Atom name for distance comparison (same group)
  --group2_atom_name CB \
# JSON list of residue-specific atom tuples, e.g. '[ ["GLY", "CA"] ]'
  --residue_specific_atoms '[["GLY", "CA"]]' \
# Include residue names in output
  --give_res_name=True \
# Include atom names in output
  --give_atom_name=True
```

```bash
# Mean distance calculation for inter-protein 

# Input GROMACS .gro structure file
mean_distance \
  --gro_file "md.gro" \
# Input GROMACS .xtc trajectory file
  --xtc_file "md.xtc" \
# GROMACS index file defining groups
  --ndx_file "index.ndx" \
# Starting frame index for analysis
  --start_frame 8000 \
# Ending frame index for analysis
  --end_frame 9000 \
# First group ID for inter-protein
  --group1_ID 10 \
# Second group ID for inter-protein
  --group2_ID 11 \
# Atom name in group (e.g., beta carbon)
  --group1_atom_name CB \
# Atom name for distance comparison (same group)
  --group2_atom_name CB \
# JSON list of residue-specific atom tuples, e.g. '[ ["GLY", "CA"] ]'
  --residue_specific_atoms '[["GLY", "CA"]]' \
# Include residue names in output
  --give_res_name=True \
# Include atom names in output
  --give_atom_name=True
```

## Parameter Notes

- `group_ID` is used **only** for intra-protein interactions analysis (single group).
- `group1_ID` and `group2_ID` are used **only** for inter-protein interactions analysis (two groups).
- `group1_atom_name` and `group2_atom_name` are **always required**:
  - For intra-protein interactions, they refer to the atom name or comma-separated atom names of interacting pairs within the **same group**.
  - For inter-protein interactions, they refer to the atom name or comma-separated atom names of interacting pairs between **different groups**.            
- You can specify the **same atom name** for both groups (e.g., `"CB"`) or **different atom names** (e.g., `"CA"` vs `"CB"`) depending on your analysis goals.
- The `residue_specific_atoms` parameter should be provided as a **Python list of tuples** (e.g., `[("GLY", "CA"), ("ALA", "CB")]`) when using the Python API, and as a **JSON-formatted string** (e.g., `'[["GLY", "CA"], ["ALA", "CB"]]'`) when using the command-line interface.
- By default, the output file is saved in the same directory as the gro_file, unless an alternative location or file path is specified via the --output_file command-line flag or the output_file="<path>" keyword argument.

## Visualising
![Example output visualization of MDInteractions showing interaction pairs over time](https://github.com/afowdar/MDInteractions/pub.png)


## License

This project is licensed under the [MIT License](LICENSE.txt).


## Citation

If you use `MDInteractions` in your research, please cite it as:

> Fowdar, A and Martin, DP. *MDInteractions: A python package for the distance-based analysis of intra- and inter-protein interactions*. GitHub, 2025. https://github.com/afowdar/MDInteractions

Or use the following BibTeX entry:

```bibtex
@misc{mdinteractions2025,
  author       = {Anjani Fowdar and Darren P Martin},
  title        = {MDInteractions: A python package for the distance-based analysis of intra- and inter-protein interactions},
  year         = {2025},
  howpublished = {\url{https://github.com/afowdar/MDInteractions}},
}
```

## Getting help

For questions, ideas, or feedback, join the conversation in the [Discussions](https://github.com/afowdar/MDInteractions/discussions) section!
Feel free to reach out: anjanifowdar@gmail.com

