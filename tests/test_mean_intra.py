import pytest
from pathlib import Path
import pandas as pd

from MDInteractions import mean_distance

# Path to test data directory
TEST_DATA_DIR = Path(__file__).parent / "test_data"

# Define a pytest fixture to set up a mean_distance analyzer for the test
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
        group_ID=11,
        group1_atom_name="CB",
        group2_atom_name="CB",
        give_res_name=True,
        give_atom_name=True,
        residue_specific_atoms=residue_atoms,
        output_file="test_output.csv"
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
        'Average_Distance'
    }

    # Check that all expected columns are present in the DataFrame
    assert expected_cols.issubset(df.columns), "Output CSV missing expected columns"

    # Define expected distance results manually for comparison
    manual_expected = {
        (5, 6, 'CB', 'CB', 'LEU', 'LEU'): 5.7,
        (5, 7, 'CA', 'CA', 'LEU', 'GLY'): 7.1,
        (5, 8, 'CB', 'CB', 'LEU', 'ALA'): 11.6,
        (6, 7, 'CA', 'CA', 'LEU', 'GLY'): 3.8,
        (6, 8, 'CB', 'CB', 'LEU', 'ALA'): 8.6,
        (7, 8, 'CA', 'CA', 'GLY', 'ALA'): 3.8,
        # Add more if available
    }

    # Convert DataFrame to dict for lookup
    result_dict = {}
    for _, row in df.iterrows():
        key = (
            row['Group1_resid'], row['Group2_resid'],
            row['Group1_atom'], row['Group2_atom'],
            row['Group1_resname'], row['Group2_resname']
        )
        result_dict[key] = row['Average_Distance']

    # Loop through all manually expected values and compare them with actual results
    for key, expected_dist in manual_expected.items():
        assert key in result_dict, f"Missing interaction {key} in output."
        actual = result_dict[key]
        assert abs(actual - expected_dist) < 0.5, (
            f"Mismatch for {key}: expected {expected_dist}, got {actual:.2f}"
        )

    # Optional cleanup
    output_path.unlink()