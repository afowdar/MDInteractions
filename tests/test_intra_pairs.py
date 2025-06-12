import pytest
from pathlib import Path
import pandas as pd

from MDInteractions import protein_interactions

# Define the path to test data
TEST_DATA_DIR = Path(__file__).parent / "test_data"

# Define a pytest fixture 
@pytest.fixture
def analyzer():
    gro_file = TEST_DATA_DIR / "adk_oplsaa.gro"
    xtc_file = TEST_DATA_DIR / "adk_oplsaa.xtc"
    ndx_file = TEST_DATA_DIR / "index_numbered.ndx"

    residue_atoms = [
        ('GLY', 'CA')
    ]

    # Return a configured instance of protein_interactions
    return protein_interactions(
        gro_file=gro_file,
        xtc_file=xtc_file,
        ndx_file=ndx_file,
        start_frame=1,
        end_frame=5,
        cutoff=7,
        group_ID=11,
        group1_atom_name="CB",
        group2_atom_name="CB",
        give_res_name=True,
        give_atom_name=True,
        residue_specific_atoms=residue_atoms,
        output_file="test_output.csv"
    )

# Define the actual test function that verifies the output of the analysis
def test_analyze_and_output(analyzer):
    analyzer.analyze()

    # Check that the output file was actually created
    output_path = analyzer.gro_file.parent / analyzer.output_file
    assert output_path.exists(), "Output CSV file was not created."

    # Load the output CSV file as a DataFrame
    df = pd.read_csv(output_path)

    # Define the set of expected column names
    expected_cols = {
        'Group1_resid', 'Group2_resid',
        'Group1_resname', 'Group2_resname',
        'Group1_atom', 'Group2_atom'
    }

    # Assert that all expected columns are present in the CSV output
    assert expected_cols.issubset(df.columns), "Output CSV missing expected columns"

    # Manually define expected residue-residue interactions
    manual_interactions = {
        (5, 6, 'CB', 'CB', 'LEU', 'LEU'),
        (6, 7, 'CA', 'CA', 'LEU', 'GLY'),
        (7, 8, 'CA', 'CA', 'GLY', 'ALA'),
        # Add more if needed
    }

    # Convert the DataFrame to a set of tuples for comparison
    result_set = set(
        zip(
            df['Group1_resid'],
            df['Group2_resid'],
            df['Group1_atom'],
            df['Group2_atom'],
            df['Group1_resname'],
            df['Group2_resname']
        )
    )

    # Verify that all manually expected interactions are found in the output
    for pair in manual_interactions:
        assert pair in result_set, f"Expected interaction {pair} not found in results."

    # Optional cleanup
    output_path.unlink()