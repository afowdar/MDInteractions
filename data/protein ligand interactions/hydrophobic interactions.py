import pandas as pd

hydrophobic = {"ALA","VAL","LEU","ILE","PRO","PHE","MET","TRP"}

df = pd.read_csv("consistent_residue_pairs.csv")

filtered = df[
    (df["Group1_resname"].isin(hydrophobic)) &
    (df["Group2_resname"].isin(hydrophobic))
]

filtered.to_csv("hydrophobic_only.csv", index=False)