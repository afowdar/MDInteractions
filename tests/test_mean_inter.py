import pytest
from pathlib import Path
import pandas as pd

from MDInteractions import mean_distance

# Path to test data directory
TEST_DATA_DIR = Path(__file__).parent / "test_data"

# This fixture sets up and returns a configured mean_distance analyzer object for testing
@pytest.fixture
def analyzer():
    gro_file = TEST_DATA_DIR / "adk_oplsaa.gro"
    xtc_file = TEST_DATA_DIR / "adk_oplsaa.xtc"
    ndx_file = TEST_DATA_DIR / "index_numbered.ndx"

    # Optionally override default atoms to use for specific residues
    residue_atoms = [
        ('GLY', 'CA')
    ]
    
    # Return a configured mean_distance object to be used in the test
    return mean_distance(
        gro_file=gro_file,
        xtc_file=xtc_file,
        ndx_file=ndx_file,
        start_frame=1,
        end_frame=5,
        group1_ID=10,
        group2_ID=11,
        group1_atom_name="CB",
        group2_atom_name="CB",
        give_res_name=True,
        give_atom_name=True,
        residue_specific_atoms=residue_atoms,
        output_file="average_distance.csv"
    )

# test function that uses the analyzer fixture
def test_mean_distance_analysis(analyzer):
    analyzer.analyze()

    # Check that the output file was actually created
    output_path = analyzer.gro_file.parent / analyzer.output_file
    assert output_path.exists(), "Output CSV file was not created."

    # Load the output CSV file as a DataFrame
    df = pd.read_csv(output_path)

    # Define the expected columns in the output CSV
    expected_cols = {
        'Group1_resid', 'Group2_resid',
        'Group1_resname', 'Group2_resname',
        'Group1_atom', 'Group2_atom',
        'Average_Distance', 'Standard_Deviation'    
    }

    # Check that all expected columns are present in the DataFrame
    assert expected_cols.issubset(df.columns), "Output CSV missing expected columns"

    # Define expected distance results manually for comparison
    manual_expected = {
    (1, 5, 'CB', 'CB', 'MET', 'LEU'): (12.2, 0.5119),
    (1, 6, 'CB', 'CB', 'MET', 'LEU'): (16.2, 0.5927),
    (1, 7, 'CA', 'CA', 'MET', 'GLY'): (18.7, 0.2362),
    (1, 8, 'CB', 'CB', 'MET', 'ALA'): (23.3, 0.4639),
    (2, 5, 'CB', 'CB', 'ARG', 'LEU'): (11.2, 0.4860),
    (2, 6, 'CB', 'CB', 'ARG', 'LEU'): (14.2, 0.3562),
    (2, 7, 'CA', 'CA', 'ARG', 'GLY'): (16.9, 0.3130),
    (2, 8, 'CB', 'CB', 'ARG', 'ALA'): (22.3, 0.3711),
    (3, 5, 'CB', 'CB', 'ILE', 'LEU'): (6.3, 0.2852),
    (3, 6, 'CB', 'CB', 'ILE', 'LEU'): (10.0, 0.1044),
    (3, 7, 'CA', 'CA', 'ILE', 'GLY'): (13.1, 0.3231),
    (3, 8, 'CB', 'CB', 'ILE', 'ALA'): (17.4, 0.2750),
    (4, 5, 'CB', 'CB', 'ILE', 'LEU'): (5.8, 0.0708),
    (4, 6, 'CB', 'CB', 'ILE', 'LEU'): (6.5, 0.0960),
    (4, 7, 'CA', 'CA', 'ILE', 'GLY'): (10.3, 0.1295),
    (4, 8, 'CB', 'CB', 'ILE', 'ALA'): (14.8, 0.1119)
    }

    # Convert DataFrame to dict for lookup
    result_dict = {}
    for _, row in df.iterrows():
        key = (
            row['Group1_resid'], row['Group2_resid'],
            row['Group1_atom'], row['Group2_atom'],
            row['Group1_resname'], row['Group2_resname']
        )
        result_dict[key] = {
        'Average_Distance': row['Average_Distance'],
        'Standard_Deviation': row['Standard_Deviation']
    } #added

    # Loop through all manually expected values and compare them with actual results
    for key, (expected_dist, expected_std) in manual_expected.items():
        assert key in result_dict, f"Missing interaction {key} in output."
        actual_dist = result_dict[key]['Average_Distance']
        actual_std = result_dict[key]['Standard_Deviation']
        # Allow a small tolerance in distance comparisons due to floating-point calculations
        assert abs(actual_dist - expected_dist) < 0.5, (
            f"Mismatch for {key}: expected distance {expected_dist}, got {actual_dist:.2f}"
        )
        assert abs(actual_std - expected_std) < 0.5, (
            f"Mismatch for {key}: expected std {expected_std}, got {actual_std:.2f}"
        )
    # Optional cleanup
    output_path.unlink()



