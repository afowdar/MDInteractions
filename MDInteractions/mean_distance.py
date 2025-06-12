from collections import defaultdict
from pathlib import Path 
import logging
from datetime import datetime

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter

logger = logging.getLogger(__name__)


class MeanDistanceAnalyzer:
    """
    Analyzes mean pairwise distances between atoms in specified groups over a trajectory.

    Supports both intra-group and inter-group analysis based on provided index groups.
    """
    def __init__(self, gro_file, xtc_file, ndx_file, start_frame, end_frame,
                 group_ID=None, group1_ID=None, group2_ID=None,
                 group1_atom_name="", group2_atom_name="",
                 give_res_name=False, give_atom_name=False,
                 residue_specific_atoms=None, output_file= None):
        """
        Initializes the MeanDistanceAnalyzer with input files and parameters.

        Args:
            gro_file (str | Path): Path to the .gro file.
            xtc_file (str | Path): Path to the .xtc file.
            ndx_file (str | Path): Path to the .ndx file.
            start_frame (int): Start frame index for analysis.
            end_frame (int): End frame index for analysis.
            group_ID (int, optional): Group ID for intra-group analysis.
            group1_ID (int, optional): Group 1 ID for inter-group analysis.
            group2_ID (int, optional): Group 2 ID for inter-group analysis.
            group1_atom_name (str): Comma-separated atom names in group 1 to consider.
            group2_atom_name (str): Comma-separated atom names in group 2 to consider.
            give_res_name (bool, optional): Whether to include residue names in the output.
            give_atom_name (bool, optional): Whether to include atom names in the output.
            residue_specific_atoms (list[tuple], optional): Custom atom names per residue
            output_file (str | Path, optional): Path to output CSV file.
        """
        
        self.gro_file = Path(gro_file)
        self.xtc_file = Path(xtc_file)
        self.ndx_file = Path(ndx_file)

        # Validate input files
        for file_path in [self.gro_file, self.xtc_file, self.ndx_file]:
            if not file_path.exists():
                logger.error(f"{file_path.name} not found!")
                raise FileNotFoundError(f"{file_path.name} not found!")

        # Determine analysis mode: intra or inter 
        self.mode = "inter" if group1_ID is not None and group2_ID is not None else "intra"
        if self.mode == "inter" and (group1_ID is None or group2_ID is None):
            raise ValueError("Both group1_ID and group2_ID must be provided for inter mode.")
        if self.mode == "intra" and group_ID is None:
            raise ValueError("group_ID must be provided for intra mode.")

        self.start_frame = start_frame
        self.end_frame = end_frame
        self.group_ids = [group1_ID, group2_ID] if self.mode == "inter" else [group_ID]
        self.give_res_name = give_res_name
        self.give_atom_name = give_atom_name

        self.group1_atom_names = [name.strip() for name in group1_atom_name.split(',')] if group1_atom_name else []
        self.group2_atom_names = [name.strip() for name in group2_atom_name.split(',')] if group2_atom_name else []

        # Build a dictionary of residue-specific atoms
        self.residue_specific_atoms = {}
        if residue_specific_atoms:
            for entry in residue_specific_atoms:
                resname = entry[0]
                atoms = entry[1:]
                self.residue_specific_atoms[resname] = atoms
                
        self.output_file = Path(output_file) if output_file is not None else None

        # These will be set in prepare()
        self.universe = None
        self.replaced_ndx_file = None
        self.atom_indices = None

    def cut_trajectory(self):
        """
        Cuts the original trajectory from start_frame to end_frame
        and writes the result to a new XTC file.
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        trimmed_xtc_path = self.xtc_file.with_name(f"{self.xtc_file.stem}_trimmed_{timestamp}.xtc")

        logger.info(f"Cutting trajectory from frame {self.start_frame} to {self.end_frame}")

        u = mda.Universe(str(self.gro_file), str(self.xtc_file))
        with XTCWriter(str(trimmed_xtc_path), n_atoms=u.atoms.n_atoms) as writer:
            for ts in u.trajectory[self.start_frame:self.end_frame + 1]:
                writer.write(u.atoms)

        logger.info(f"Trimmed trajectory saved to {trimmed_xtc_path}")
        return trimmed_xtc_path

    def replace_group_names(self):
        """
        Creates a new .ndx file with numeric group labels.
        """
        logger.debug("Replacing group names in .ndx file")
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

        base_name = self.ndx_file.stem
        output_file = self.ndx_file.parent / f"{base_name}_numbered.ndx"

        with output_file.open('w') as f:
            f.writelines(output_lines)

        return output_file

    def get_group_atom_indices(self):
        """
        Parses the numeric-labeled NDX file to extract atom indices per group.
        """
        logger.debug("Extracting atom indices from replaced .ndx file")
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
        Prepares the trimmed trajectory and loads the necessary atom indices.
        """
        trimmed_xtc = self.cut_trajectory()
        self.universe = mda.Universe(str(self.gro_file), str(trimmed_xtc))
        self.replaced_ndx_file = self.replace_group_names()
        self.atom_indices = self.get_group_atom_indices()

    def compare_residues(self, res1, res2, pairwise_distance):
        """
        Records distances between specified atom pairs.

        Supports optional residue-specific atom matching or falls back to user-defined atom names.
        """
        override_atoms = None
        if res1.resname in self.residue_specific_atoms:
            override_atoms = self.residue_specific_atoms[res1.resname]
        elif res2.resname in self.residue_specific_atoms:
            override_atoms = self.residue_specific_atoms[res2.resname]

        if override_atoms:
            atom1_names = atom2_names = override_atoms
        else:
            atom1_names = self.group1_atom_names
            atom2_names = self.group2_atom_names 
        
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
                distance = np.linalg.norm(atom1.position - atom2.position)
                key = ((res1.resid, res1.resname), (res2.resid, res2.resname), atom1_name, atom2_name)
                pairwise_distance[key].append(distance)

    def process_frame(self, ts):
        """
        Processes a single trajectory frame and returns pairwise distances.
        """
        pairwise_distance = defaultdict(list)

        try:
            if self.mode == 'intra':
                atoms = self.universe.select_atoms(
                    f"bynum {' '.join(map(str, self.atom_indices[self.group_ids[0]]))}")
                residue_list = list(atoms.residues)
                for i, res1 in enumerate(residue_list):
                    for res2 in residue_list[i + 1:]:
                        self.compare_residues(res1, res2, pairwise_distance)

            elif self.mode == 'inter':
                atoms1 = self.universe.select_atoms(
                    f"bynum {' '.join(map(str, self.atom_indices[self.group_ids[0]]))}")
                atoms2 = self.universe.select_atoms(
                    f"bynum {' '.join(map(str, self.atom_indices[self.group_ids[1]]))}")
                for res1 in atoms1.residues:
                    for res2 in atoms2.residues:
                        self.compare_residues(res1, res2, pairwise_distance)

        except Exception as e:
            logger.warning(f"Error during frame {ts.frame} processing: {e}")

        return pairwise_distance

    def process_frame_by_index(self, frame_index):
        """
        Allows frame processing externally (e.g., for multiprocessing).
        """
        ts = self.universe.trajectory[frame_index]
        return self.process_frame(ts)
    
    def aggregate_distances(self, all_distances):
        """
        Aggregates all pairwise distances across frames and writes output.
        """
        combined_counts = defaultdict(list)
        for pairwise in all_distances:
            for key, dists in pairwise.items():
                combined_counts[key].extend(dists)

        x, y, x_name, y_name, x_atoms, y_atoms, avg_dist = [], [], [], [], [], [], []

        for key, distances in combined_counts.items():
            (resid1, resname1), (resid2, resname2), atom1_name, atom2_name = key
            avg_distance = np.mean(distances)
            x.append(resid1)
            y.append(resid2)
            avg_dist.append(avg_distance)

            if self.give_res_name:
                x_name.append(resname1)
                y_name.append(resname2)
            if self.give_atom_name:
                x_atoms.append(atom1_name)
                y_atoms.append(atom2_name)

        # Construct output DataFrame
        data = {
            'Group1_resid': x,
            'Group2_resid': y,
            'Average_Distance': avg_dist
        }
        if self.give_res_name:
            data['Group1_resname'] = x_name
            data['Group2_resname'] = y_name
        if self.give_atom_name:
            data['Group1_atom'] = x_atoms
            data['Group2_atom'] = y_atoms

        df = pd.DataFrame(data)
        filename = self.output_file or (
            "average_distance.csv" if self.mode == "inter"
            else f"average_distance_group{self.group_ids[0]}.csv"
        )
        output_path = self.gro_file.parent / filename
        df.to_csv(output_path, index=False)
        logger.info(f"Analysis complete. Results saved to {output_path}")

    def analyze(self):
        """
        Runs the full analysis pipeline: prepare trajectory, process all frames, and aggregate.
        """
        self.prepare()
        all_distances = []

        for ts in self.universe.trajectory:
            frame_distances = self.process_frame(ts)
            all_distances.append(frame_distances)

        self.aggregate_distances(all_distances)


def mean_distance(*args, **kwargs):
    """
    A convenience wrapper function to instantiate MeanDistanceAnalyzer.
    """
    return MeanDistanceAnalyzer(*args, **kwargs)


__all__ = ['mean_distance']
