import pandas as pd

# Load Sample Dataset (Case Group)
sample_df = pd.read_csv("combined_RNA-seq_data.csv", index_col=0)

# Load Control Dataset
control_df = pd.read_csv("filtered_control_data.csv", index_col=0)

# Debug: Print first few rows and column names to check structure
print("Control DataFrame Columns:", control_df.columns)

# Check if 'Description' exists, then rename it to 'gene_name'
if 'Description' in control_df.columns:
    control_df.rename(columns={'Description': 'gene_name'}, inplace=True)

# Debug: Print again to confirm renaming
print("Renamed Control DataFrame Columns:", control_df.columns)

# Ensure 'gene_name' column exists before merging
if 'gene_name' in control_df.columns:
    merged_df = sample_df.merge(control_df[['gene_name']], left_index=True, right_on='gene_name')
    merged_df.set_index("gene_name", inplace=True)
    print(f"Number of overlapping genes: {merged_df.shape[0]}")
else:
    print("ERROR: 'gene_name' column not found in control data. Check dataset format.")
