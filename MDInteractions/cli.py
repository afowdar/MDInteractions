import os
import ast
import argparse
import logging
from pathlib import Path
from collections import defaultdict
from functools import partial
import concurrent.futures

import numpy as np
import MDAnalysis as mda

from MDInteractions.protein_interactions import protein_interactions

# Set up logger for this module
logger = logging.getLogger(__name__)


def worker_process_frame(frame_index, gro_file, xtc_file, group_ids, atom_indices,
                         group1_atom_names, group2_atom_names, mode, residue_specific_atoms,
                         cutoff_min, cutoff_max):
    """
    Process a single frame of the trajectory to identify residue-residue interactions
    within a specified cutoff range.

    Parameters:
        frame_index (int): The index of the trajectory frame to analyze.
        gro_file (Path): Path to the input .gro structure file.
        xtc_file (Path): Path to the input .xtc trajectory file.
        group_ids (list): Group ID(s) specifying which index groups to analyze.
        atom_indices (dict): Dictionary mapping group IDs to atom indices.
        group1_atom_names (list): Atom name or comma-separated atom names for first atom in the interacting pair.
        group2_atom_names (list): Atom name or comma-separated atom names for second atom in the interacting pair.
        mode (str): Either 'inter' or 'intra' specifying interaction mode.
        residue_specific_atoms (list) or None): Optional override of atom names by residue name.
        cutoff_min (float): Minimum distance for interaction.
        cutoff_max (float): Maximum distance for interaction.
    """

    # Load the universe and select the specified frame
    universe = mda.Universe(str(gro_file), str(xtc_file))
    universe.trajectory[frame_index]

    # Initialize results container
    pairs_found = defaultdict(int)

    def compare_residues(res1, res2):
        """
        Compare a pair of residues to identify atom-atom contacts within the distance cutoff.
        """
        # Determine atom names to use (either overridden per-residue or general group lists)
        override_atoms = None
        if residue_specific_atoms:
            if res1.resname in residue_specific_atoms:
                override_atoms = residue_specific_atoms[res1.resname]
            elif res2.resname in residue_specific_atoms:
                override_atoms = residue_specific_atoms[res2.resname]

        atom1_names, atom2_names = (override_atoms, override_atoms) if override_atoms else (group1_atom_names, group2_atom_names)

        # Iterate through atom pairs and compute distances
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
                if cutoff_min < dist <= cutoff_max:
                    key = ((res1.resid, res1.resname), (res2.resid, res2.resname), atom1_name, atom2_name)
                    pairs_found[key] += 1

    # Process residues based on mode (intra-group or inter-group)
    try:
        if mode == 'intra':
            atoms = universe.select_atoms(f"bynum {' '.join(map(str, atom_indices[group_ids[0]]))}")
            residues = list(atoms.residues)
            for i, res1 in enumerate(residues):
                for res2 in residues[i+1:]:
                    compare_residues(res1, res2)

        elif mode == 'inter':
            atoms1 = universe.select_atoms(f"bynum {' '.join(map(str, atom_indices[group_ids[0]]))}")
            atoms2 = universe.select_atoms(f"bynum {' '.join(map(str, atom_indices[group_ids[1]]))}")
            for res1 in atoms1.residues:
                for res2 in atoms2.residues:
                    compare_residues(res1, res2)
    except Exception as e:
        logger.warning(f"Error in frame {frame_index}: {e}")
    return pairs_found

def parse_args():
    """
    Parse command-line arguments for the interaction analysis script.
    """
    parser = argparse.ArgumentParser(description="Analyze protein residue-residue interactions from MD trajectories.")

    # Required input files
    parser.add_argument("--gro_file", type=Path, required=True, help="Path to the .gro file")
    parser.add_argument("--xtc_file", type=Path, required=True, help="Path to the .xtc file")
    parser.add_argument("--ndx_file", type=Path, required=True, help="Path to the .ndx file")
    
    # Trajectory frame range
    parser.add_argument("--start_frame", type=int, required=True, help="Start frame for analysis")
    parser.add_argument("--end_frame", type=int, required=True, help="End frame for analysis")

    # Group selection
    parser.add_argument("--group_ID", type=int, help="Group ID for intra mode (single group).")
    parser.add_argument("--group1_ID", type=int, help="Group 1 ID for inter mode.")
    parser.add_argument("--group2_ID", type=int, help="Group 2 ID for inter mode.")

    # Atom selection
    parser.add_argument("--group1_atom_name", type=str, default="", help="Atom name or comma-separated atom names for first atom in the interacting pair.")
    parser.add_argument("--group2_atom_name", type=str, default="", help="Atom name or comma-separated atom names for second atom in the interacting pair.")

    # Distance cutoff specification
    parser.add_argument("--cutoff", type=str, required=True,
                        help="Cutoff distance or range, e.g. '5.0' or '[3.0, 5.0]'")

    # Output options
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

    # Optional atom-specific override
    parser.add_argument(
        "--residue_specific_atoms",
        type=str,
        default=None,
        help="Python-style list, e.g. '[ [\"GLY\", \"CA\"] ]'"
    )

    # Output file path
    parser.add_argument(
        "--output_file",
        type=Path,
        default=None,
        help="Optional output CSV filename"
    )

    # Logging and parallelism
    parser.add_argument(
        "--log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level"
    )

    parser.add_argument(
        "--num_workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: all cores or SLURM_CPUS_PER_TASK)"
    )
    return parser.parse_args()


def main():
    """
    Main function to perform residue-residue interaction analysis over a molecular dynamics trajectory.
    
    Steps:
        - Parses command-line arguments.
        - Configures logging.
        - Safely interprets cutoff and residue-specific atom settings.
        - Initializes the `protein_interactions` analyzer.
        - Prepares the MDAnalysis universe and index groups.
        - Launches parallel processing over trajectory frames.
        - Aggregates and optionally saves the interaction results.
    """
    args = parse_args()

    # Set up logging format and level
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)

    # Safely evaluate cutoff input (single float or list of two floats)
    try:
        cutoff_val = ast.literal_eval(args.cutoff)
        if isinstance(cutoff_val, (int, float)):
            cutoff_min = 0.0
            cutoff_max = float(cutoff_val)
        elif isinstance(cutoff_val, (list, tuple)) and len(cutoff_val) == 2:
            cutoff_min = min(cutoff_val)
            cutoff_max = max(cutoff_val)
        else:
            raise ValueError("Cutoff must be a number or list/tuple of two numbers.")
    except Exception as e:
        logger.error(f"Invalid cutoff argument: {e}")
        return

    # Safely evaluate residue-specific atom overrides if provided
    residue_specific_atoms = None
    if args.residue_specific_atoms:
        try:
            parsed = ast.literal_eval(args.residue_specific_atoms)
            if isinstance(parsed, list) and all(isinstance(pair, (list, tuple)) and len(pair) == 2 for pair in parsed):
                residue_specific_atoms = {}
                for res, atom in parsed:
                    residue_specific_atoms.setdefault(res, []).append(atom)
            else:
                raise ValueError("residue_specific_atoms must be list of pairs")
        except Exception as e:
            logger.error(f"Failed to parse residue_specific_atoms: {e}")
            return

    # Initialize the protein interaction analyzer
    analyzer = protein_interactions(
        gro_file=args.gro_file,
        xtc_file=args.xtc_file,
        ndx_file=args.ndx_file,
        start_frame=args.start_frame,
        end_frame=args.end_frame,
        cutoff=(cutoff_min, cutoff_max),
        group_ID=args.group_ID,
        group1_ID=args.group1_ID,
        group2_ID=args.group2_ID,
        group1_atom_name=args.group1_atom_name,
        group2_atom_name=args.group2_atom_name,
        give_res_name=args.give_res_name,
        give_atom_name=args.give_atom_name,
        residue_specific_atoms=residue_specific_atoms,
        output_file=args.output_file,
    )

    logger.info("Preparing analyzer (loading trajectory and index groups)...")
    analyzer.prepare()

    # Make sure resid_to_resname attribute exists in analyzer for output
    if args.give_res_name and not hasattr(analyzer, "resid_to_resname"):
        analyzer.resid_to_resname = {}
        for res in analyzer.universe.residues:
            analyzer.resid_to_resname[res.resid] = res.resname

    # Determine number of workers (parallel processes) to use
    # First try user-specified, then SLURM environment variable, finally system CPU coun
    num_workers = args.num_workers or int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count()))
    logger.info(f"Using {num_workers} parallel workers")

    frame_indices = [ts.frame for ts in analyzer.universe.trajectory]

    worker_func = partial(
        worker_process_frame,
        gro_file=analyzer.gro_file,
        xtc_file=analyzer.universe.trajectory.filename,
        group_ids=analyzer.group_ids,
        atom_indices=analyzer.atom_indices,
        group1_atom_names=analyzer.group1_atom_names,
        group2_atom_names=analyzer.group2_atom_names,
        mode=analyzer.mode,
        residue_specific_atoms=residue_specific_atoms,
        cutoff_min=cutoff_min,
        cutoff_max=cutoff_max,
    )

    # Use multiprocessing to process each frame in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Execute worker_func on each frame index and collect results
        all_results = list(executor.map(worker_func, frame_indices))

    # Aggregate results from all frames into a final output
    analyzer.aggregate_results(all_results)

if __name__ == "__main__":
    main()
