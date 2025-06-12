from collections import defaultdict
from pathlib import Path 
import logging
from datetime import datetime

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter

logger = logging.getLogger(__name__)


class ProteinInteractionAnalyzer:
    """
    Analyzes protein residue-residue interactions from molecular dynamics trajectories using MDAnalysis.
    Supports both intra-group (within the same selection) and inter-group (between two selections) residue
    interaction detection with support for atom-level specificity and residue-specific atom mappings.
    """

    def __init__(self, gro_file, xtc_file, ndx_file, start_frame: int, end_frame: int,
                 cutoff: float | list[float] | tuple[float, float],
                 group_ID: int = None, group1_ID: int = None, group2_ID: int = None,
                 group1_atom_name: str = "", group2_atom_name: str = "",
                 give_res_name: bool = False, give_atom_name: bool = False,
                 residue_specific_atoms=None, output_file= None):
        """
        Initializes the ProteinInteractionAnalyzer with necessary files, parameters, and MDAnalysis universe.
        Validates the existence of input files and determines whether intra- or inter-group analysis will be performed.

        Args:
            gro_file (str | Path): Path to the .gro file.
            xtc_file (str | Path): Path to the .xtc file.
            ndx_file (str | Path): Path to the .ndx file.
            start_frame (int): Starting frame index.
            end_frame (int): Ending frame index.
            cutoff (float | list[float] | tuple[float, float]): Distance cutoff or range.
            group_ID (int, optional): Group ID for intra analysis.
            group1_ID (int, optional): Group ID for group 1 in inter analysis.
            group2_ID (int, optional): Group ID for group 2 in inter analysis.
            group1_atom_name (str): Atom name or comma-separated atom names for first atom in the interacting pair.
            group2_atom_name (str): Atom name or comma-separated atom names for second atom in the interacting pair.
            give_res_name (bool, optional): Whether to output residue names.
            give_atom_name (bool, optional): Whether to output atom names.
            residue_specific_atoms (list[tuple], optional): Custom atom names per residue.
            output_file (str | Path, optional): Filename to save the output.
        """
        
        # Convert file paths to Path objects
        self.gro_file = Path(gro_file)
        self.xtc_file = Path(xtc_file)
        self.ndx_file = Path(ndx_file)

        # Validate input files
        for file_path in [self.gro_file, self.xtc_file, self.ndx_file]:
            if not file_path.exists():
                logger.error(f"{file_path.name} not found!")
                raise FileNotFoundError(f"{file_path.name} not found!")
            
        # Determine mode of analysis (intra or inter)
        self.mode = "inter" if group1_ID is not None and group2_ID is not None else "intra"
        if self.mode == "inter" and (group1_ID is None or group2_ID is None):
            raise ValueError("Both group1_ID and group2_ID must be provided for inter mode.")
        if self.mode == "intra" and group_ID is None:
            raise ValueError("group_ID must be provided for intra mode.")

        # Frame range for trimming trajectory
        self.start_frame = start_frame
        self.end_frame = end_frame
        
        self.group_ids = [group1_ID, group2_ID] if self.mode == "inter" else [group_ID]
        
        # Output format options
        self.give_res_name = give_res_name
        self.give_atom_name = give_atom_name

        # Atom names per group
        self.group1_atom_names = [name.strip() for name in group1_atom_name.split(',')] if group1_atom_name else []
        self.group2_atom_names = [name.strip() for name in group2_atom_name.split(',')] if group2_atom_name else []

        # Process distance cutoff
        if isinstance(cutoff, (int, float)):
            self.cutoff_min = 0.0
            self.cutoff_max = cutoff
        elif isinstance(cutoff, (list, tuple)) and len(cutoff) == 2:
            self.cutoff_min = min(cutoff)
            self.cutoff_max = max(cutoff)
        else:
            raise ValueError("cutoff must be a number or a list/tuple of two numbers")

        # Custom residue-atom mapping if provided
        self.residue_specific_atoms = {}
        if residue_specific_atoms:
            for entry in residue_specific_atoms:
                resname = entry[0]
                atoms = entry[1:]
                self.residue_specific_atoms[resname] = atoms

        # Optional output file path
        self.output_file = Path(output_file) if output_file is not None else None

        # Internal variables
        self.universe = None
        self.replaced_ndx_file = None
        self.atom_indices = None

    def cut_trajectory(self):
        """
        Trim trajectory to the specified frame range and save as new .xtc file.
        Helps to reduce memory and speed up processing.
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        trimmed_xtc_path = self.xtc_file.with_name(f"{self.xtc_file.stem}_trimmed_{timestamp}.xtc")

        logger.info(f"Cutting trajectory from frame {self.start_frame} to {self.end_frame}")
        u = mda.Universe(self.gro_file, self.xtc_file)
        
        # Write only the desired frames
        with XTCWriter(str(trimmed_xtc_path), n_atoms=u.atoms.n_atoms) as writer:
            for ts in u.trajectory[self.start_frame:self.end_frame + 1]:
                writer.write(u.atoms)

        logger.info(f"Trimmed trajectory written to {trimmed_xtc_path}")
        return trimmed_xtc_path

    def replace_group_names(self):
        """
        Replace named groups in the .ndx file with numeric identifiers.
        """
        group_counter = 0
        output_lines = []
        with self.ndx_file.open('r') as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith('[') and stripped.endswith(']'):
                    output_lines.append(f"[ {group_counter} ]\n")
                    group_counter += 1
                else:
                    output_lines.append(line)

        output_file = self.ndx_file.with_name(self.ndx_file.stem + "_numbered.ndx")
        with output_file.open('w') as f:
            f.writelines(output_lines)
        return output_file

    def get_group_atom_indices(self):
        """
        Parses the numeric .ndx file and extracts atom indices for each group.
        """
        group_atoms = {}
        with self.replaced_ndx_file.open('r') as f:
            lines = f.readlines()

        reading_group = None
        for line in lines:
            if line.startswith('[') and line.endswith(']\n'):
                try:
                    group_number = int(line.strip().strip('[]'))
                    reading_group = group_number
                    group_atoms[reading_group] = []
                except ValueError:
                    reading_group = None
            elif reading_group is not None:
                group_atoms[reading_group].extend(map(int, line.split()))

        return group_atoms

    def prepare(self):
        """
        Prepares the trimmed trajectory and loads atom group indices from the .ndx file.
        This step is required before running the analysis.
        """
        logger.info("Preparing analysis: trimming trajectory and processing index file.")
        trimmed_xtc = self.cut_trajectory()
        self.universe = mda.Universe(self.gro_file, trimmed_xtc)
        self.replaced_ndx_file = self.replace_group_names()
        self.atom_indices = self.get_group_atom_indices()

    def compare_residues(self, res1, res2, consistent_pairs):
        """
        Determine whether specific atoms of two residues are within cutoff range.
        Records residue pair if condition is met.

        Supports optional residue-specific atom matching or falls back to user-defined atom names.
        """
        override_atoms = None
        if res1.resname in self.residue_specific_atoms:
            override_atoms = self.residue_specific_atoms[res1.resname]
        elif res2.resname in self.residue_specific_atoms:
            override_atoms = self.residue_specific_atoms[res2.resname]

        atom1_names, atom2_names = (
            (override_atoms, override_atoms) if override_atoms else
            (self.group1_atom_names, self.group2_atom_names)
        )

        for atom1_name in atom1_names:
            atom1_sel = res1.atoms.select_atoms(f"name {atom1_name}")
            if len(atom1_sel) == 0:
                continue
            atom1 = atom1_sel[0]

            for atom2_name in atom2_names:
                atom2_sel = res2.atoms.select_atoms(f"name {atom2_name}")
                if len(atom2_sel) == 0:
                    continue
                atom2 = atom2_sel[0]

                # Compute Euclidean distance between atoms
                distance = np.linalg.norm(atom1.position - atom2.position)
                if self.cutoff_min < distance <= self.cutoff_max:
                    key = ((res1.resid, res1.resname), (res2.resid, res2.resname), atom1_name, atom2_name)
                    consistent_pairs[key] += 1

    def process_frame(self, ts):
        """
        Analyze a single frame in the trajectory and return the set of interacting residue pairs.
        """
        pairs_found = defaultdict(int)
        try:
            if self.mode == 'intra':
                atoms = self.universe.select_atoms(f"bynum {' '.join(map(str, self.atom_indices[self.group_ids[0]]))}")
                residue_list = list(atoms.residues)
                for i, res1 in enumerate(residue_list):
                    for j in range(i + 1, len(residue_list)):
                        res2 = residue_list[j]
                        self.compare_residues(res1, res2, pairs_found)

            elif self.mode == 'inter':
                atoms1 = self.universe.select_atoms(f"bynum {' '.join(map(str, self.atom_indices[self.group_ids[0]]))}")
                atoms2 = self.universe.select_atoms(f"bynum {' '.join(map(str, self.atom_indices[self.group_ids[1]]))}")
                for res1 in atoms1.residues:
                    for res2 in atoms2.residues:
                        self.compare_residues(res1, res2, pairs_found)
        except Exception as e:
            logger.warning(f"Error during frame {ts.frame} processing: {e}")
        return pairs_found

    def process_frame_by_index(self, frame_index):
        """
        Allows parallel frame processing by frame index.
        Used for multiprocessing from CLI or other scripts.
        """
        ts = self.universe.trajectory[frame_index]
        return self.process_frame(ts)

    def aggregate_results(self, all_pairs_counts):
        """
        Aggregates interaction results across all frames and writes consistent interactions to a CSV file.
        """
        total_frames = len(all_pairs_counts)
        combined_counts = defaultdict(int)
        for frame_dict in all_pairs_counts:
            for key, count in frame_dict.items():
                combined_counts[key] += count

        # Keep only interactions that occur in every frame
        consistent_pairs = {k: v for k, v in combined_counts.items() if v == total_frames}

        # Prepare DataFrame
        x, y, x_name, y_name, x_atoms, y_atoms = [], [], [], [], [], []
        for ((resid1, resname1), (resid2, resname2), atom1_name, atom2_name), count in consistent_pairs.items():
            x.append(resid1)
            y.append(resid2)
            if self.give_res_name:
                x_name.append(resname1)
                y_name.append(resname2)
            if self.give_atom_name:
                x_atoms.append(atom1_name)
                y_atoms.append(atom2_name)

        data = {'Group1_resid': x, 'Group2_resid': y}
        if self.give_res_name:
            data['Group1_resname'] = x_name
            data['Group2_resname'] = y_name
        if self.give_atom_name:
            data['Group1_atom'] = x_atoms
            data['Group2_atom'] = y_atoms

        df = pd.DataFrame(data)
        filename = self.output_file or (
            "consistent_residue_pairs.csv" if self.mode == "inter"
            else f"consistent_residue_pairs_group{self.group_ids[0]}.csv"
        )
        output_path = self.gro_file.parent / filename
        df.to_csv(output_path, index=False)
        logger.info(f"Analysis complete. Results saved to {output_path}")
    
    def analyze(self):
        """
        Full analysis workflow: prepares input, processes each frame, and aggregates consistent interactions.
        """
        self.prepare()
        all_pairs_counts = []
        for ts in self.universe.trajectory:
            pairs = self.process_frame(ts)
            all_pairs_counts.append(pairs)
        self.aggregate_results(all_pairs_counts)


def protein_interactions(*args, **kwargs):
    """
    A convenience wrapper function to instantiate ProteinInteractionAnalyzer.
    """
    return ProteinInteractionAnalyzer(*args, **kwargs)


__all__ = ['protein_interactions']
