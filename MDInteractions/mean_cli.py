import os
import ast
import argparse
import logging
from pathlib import Path
from collections import defaultdict
import concurrent.futures
import logging

import numpy as np
import MDAnalysis as mda

from MDInteractions.mean_distance import mean_distance

# Set up a module-level logger
logger = logging.getLogger(__name__)

def worker_process_frame(frame_index, gro_file, xtc_file, group_ids, atom_indices,
                         group1_atom_name, group2_atom_name, mode, residue_specific_atoms):
    """
    Worker function to process a single frame of a molecular dynamics trajectory.
    
    This function is designed to be executed in parallel using multiprocessing. It reconstructs
    the MDAnalysis Universe for the specified frame and computes pairwise distances between
    atoms of interest based on provided selection criteria.

    Args:
        frame_index (int): Index of the trajectory frame to process.
        gro_file (str or Path): Path to the GRO structure file.
        xtc_file (str or Path): Path to the XTC trajectory file.
        group_ids (list): List of group IDs for residue selection (1 or 2 groups depending on mode).
        atom_indices (dict): Mapping from group IDs to lists of atom indices.
        group1_atom_name (list): List of atom names for group 1.
        group2_atom_name (list): List of atom names for group 2.
        mode (str): Analysis mode, either 'intra' (within group) or 'inter' (between groups).
        residue_specific_atoms (list or None): Optional override of atom names per residue.
    """
    # Create Universe and set the frame to the one being processed
    universe = mda.Universe(str(gro_file), str(xtc_file))
    universe.trajectory[frame_index]
    pairwise_distance = defaultdict(list)

    def compare_residues(res1, res2):
        # Determine atom names to use based on residue-specific overrides
        override_atoms = None
        if residue_specific_atoms:
            if res1.resname in residue_specific_atoms:
                override_atoms = residue_specific_atoms[res1.resname]
            elif res2.resname in residue_specific_atoms:
                override_atoms = residue_specific_atoms[res2.resname]

        # Use override if available, otherwise use group-specific names
        if override_atoms:
            atom1_names = atom2_names = override_atoms
        else:
            atom1_names = group1_atom_name
            atom2_names = group2_atom_name

        # Calculate distances between all specified atom pairs
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
                dist = np.linalg.norm(atom1.position - atom2.position)
                key = ((res1.resid, res1.resname), (res2.resid, res2.resname), atom1_name, atom2_name)
                pairwise_distance[key].append(dist)

    try:
        # INTRA-mode: distances within the same group
        if mode == 'intra':
            atoms = universe.select_atoms(f"bynum {' '.join(map(str, atom_indices[group_ids[0]]))}")
            residues = list(atoms.residues)
            for i, res1 in enumerate(residues):
                for res2 in residues[i+1:]:
                    compare_residues(res1, res2)

        # INTER-mode: distances between two groups
        elif mode == 'inter':
            atoms1 = universe.select_atoms(f"bynum {' '.join(map(str, atom_indices[group_ids[0]]))}")
            atoms2 = universe.select_atoms(f"bynum {' '.join(map(str, atom_indices[group_ids[1]]))}")
            for res1 in atoms1.residues:
                for res2 in atoms2.residues:
                    compare_residues(res1, res2)
    except Exception as e:
        logger.warning(f"Error in frame {frame_index}: {e}")

    return pairwise_distance


def parse_args():
    """
    Parses command-line arguments for the MDInteractions analysis tool.
    """
    parser = argparse.ArgumentParser(description="Analyze average residue-residue distances from MD trajectories.")
    
    # Input files
    parser.add_argument("--gro_file", type=Path, required=True, help="Path to the .gro file")
    parser.add_argument("--xtc_file", type=Path, required=True, help="Path to the .xtc file")
    parser.add_argument("--ndx_file", type=Path, required=True, help="Path to the .ndx file")
    
    # Frame range
    parser.add_argument("--start_frame", type=int, required=True, help="Start frame for analysis")
    parser.add_argument("--end_frame", type=int, required=True, help="End frame for analysis")

    # Group IDs and atom selection
    parser.add_argument("--group_ID", type=int, help="Group ID for intra mode (single group).")
    parser.add_argument("--group1_ID", type=int, help="Group 1 ID for inter mode.")
    parser.add_argument("--group2_ID", type=int, help="Group 2 ID for inter mode.")
    parser.add_argument("--group1_atom_name", type=str, default="", help="Atom name or comma-separated atom names for first atom in the interacting pair.")
    parser.add_argument("--group2_atom_name", type=str, default="", help="Atom name or comma-separated atom names for second atom in the interacting pair.")

    # Output configuration
    parser.add_argument(
        "--give_res_name",
        type=lambda x: x.lower() in ["true", "1", "yes"],
        default=False,
        help="Include residue names in output."
    )
    parser.add_argument(
        "--give_atom_name",
        type=lambda x: x.lower() in ["true", "1", "yes"],
        default=False,
        help="Include atom names in output."
    )

    # Atom name overrides per residue
    parser.add_argument(
        "--residue_specific_atoms",
        type=str,
        default=None,
        help="Python-style list of [residue, atom] pairs, e.g., '[ [\"GLY\", \"CA\"] ]'."
    )

    # Optional output file
    parser.add_argument(
        "--output_file",
        type=Path,
        default=None,
        help="Optional output CSV filename"
    )

    # Logging and parallelism
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO)"
    )
    parser.add_argument(
        "--num_workers",
        type=int,
        default=None,
        help="Number of parallel worker processes (default: all cores)"
    )
    return parser.parse_args()


def main():
    """
    Main entry point for the mean distance analysis.

    Parses arguments, configures logging, initializes the analysis class, executes
    parallel frame-wise mean distance calculations, and aggregates the results.
    """
    args = parse_args()

    # Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)
    logger.setLevel(getattr(logging, args.log_level.upper()))

    # Parse residue-specific atoms if provided
    residue_specific_atoms = None
    if args.residue_specific_atoms:
        try:
            parsed = ast.literal_eval(args.residue_specific_atoms)
            if isinstance(parsed, list) and all(isinstance(pair, (list, tuple)) and len(pair) == 2 for pair in parsed):
                residue_specific_atoms = {}
                for res, atom in parsed:
                    if res not in residue_specific_atoms:
                        residue_specific_atoms[res] = []
                    residue_specific_atoms[res].append(atom)
            else:
                raise ValueError("residue_specific_atoms must be a list of [res, atom] pairs.")
        except Exception as e:
            logger.error(f"Failed to parse residue_specific_atoms: {e}")
            raise

    # Determine number of workers
    if args.num_workers is not None:
        num_workers = args.num_workers
    else:
        num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count()))

    # Initialize analysis object
    analyzer = mean_distance(
        gro_file=args.gro_file,
        xtc_file=args.xtc_file,
        ndx_file=args.ndx_file,
        start_frame=args.start_frame,
        end_frame=args.end_frame,
        group_ID=args.group_ID,
        group1_ID=args.group1_ID,
        group2_ID=args.group2_ID,
        group1_atom_name=args.group1_atom_name,
        group2_atom_name=args.group2_atom_name,
        residue_specific_atoms=residue_specific_atoms,
        give_res_name=args.give_res_name,
        give_atom_name=args.give_atom_name,
        output_file=args.output_file
    )

    # Prepares the trimmed trajectory and loads the necessary atom indices.
    analyzer.prepare()

    logger.info(f"Starting parallel frame analysis using {num_workers} workers...")

    # Get frame indices from the loaded trajectory
    frame_indices = [ts.frame for ts in analyzer.universe.trajectory]

    # Prepare a worker function with fixed arguments
    from functools import partial
    worker_func = partial(
        worker_process_frame,
        gro_file=str(args.gro_file),
        xtc_file=str(analyzer.universe.trajectory.filename),  # trimmed xtc path
        group_ids=analyzer.group_ids,
        atom_indices=analyzer.atom_indices,
        group1_atom_name=analyzer.group1_atom_names,
        group2_atom_name=analyzer.group2_atom_names,
        mode=analyzer.mode,
        residue_specific_atoms=residue_specific_atoms
    )
    
    # Use multiprocessing to process each frame in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Execute worker_func on each frame index and collect results
        results = list(executor.map(worker_func, frame_indices))

    # Aggregate results from all frames into a final output
    analyzer.aggregate_distances(results)

if __name__ == "__main__":
    main()
