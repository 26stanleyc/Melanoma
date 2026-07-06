import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

# Load Sample Dataset (Case Group)
sample_df = pd.read_csv("/Users/stanleychen/git/Melanoma/final_data/RNA_cancer_related_genes.csv", index_col=0)
print("Sample DataFrame Columns:", sample_df.columns)
print(sample_df.head())

# Load Control Dataset
control_df = pd.read_csv("/Users/stanleychen/git/Melanoma/GSEA/filtered_control_data.csv", sep="\t")
print("Control DataFrame Columns:", control_df.columns)
print(control_df.head())

# The key issue: Control dataset uses Ensembl IDs as 'Name' and gene symbols as 'Description'
# We need to use the 'Description' column from control_df to match with sample_df's index (gene symbols)

# Create a copy of control_df with gene symbols as the index
# First check if any Description values are duplicated (this would cause issues)
duplicated_genes = control_df['Description'].duplicated()
if duplicated_genes.any():
    print(f"Warning: {duplicated_genes.sum()} duplicated gene symbols found in control dataset.")
    print("Examples:", control_df.loc[duplicated_genes, 'Description'].head())
    # Option 1: Keep first occurrence of each gene
    control_df = control_df.drop_duplicates(subset='Description', keep='first')
    print(f"After removing duplicates, control dataset has {control_df.shape[0]} rows")

# Set Description (gene symbols) as the index for control data
control_df = control_df.set_index('Description')

# Drop the 'Name' column (Ensembl IDs) since we're using gene symbols
if 'Name' in control_df.columns:
    control_df = control_df.drop('Name', axis=1)

# Now find overlapping genes (common indices)
common_genes = sample_df.index.intersection(control_df.index)
print(f"Number of overlapping genes: {len(common_genes)}")

# Filter both dataframes to include only common genes
sample_filtered = sample_df.loc[common_genes]
control_filtered = control_df.loc[common_genes]

# Concatenate the dataframes horizontally (axis=1)
combined_df = pd.concat([sample_filtered, control_filtered], axis=1)

print(f"Combined dataframe shape: {combined_df.shape}")
print("Combined dataframe sample:")
print(combined_df.head())

# Optionally save the combined dataframe
combined_df.to_csv("/Users/stanleychen/git/Melanoma/combined_case_control.csv")



# # Extract sample and control sample names
# case_samples = sample_df.columns.tolist()  # All sample columns from case dataset
# control_samples = control_df.columns.tolist()[1:]  # Exclude 'gene_name' column

# # Compute mean expression for cases and controls
# case_means = merged_df[case_samples].mean(axis=1)
# control_means = control_df.set_index("gene_name").loc[merged_df.index, control_samples].mean(axis=1)

# # Compute Log-Fold Change (LFC)
# log_fc = np.log2(case_means + 1) - np.log2(control_means + 1)

# # Perform t-test to get p-values
# ttest_pvalues = ttest_ind(
#     merged_df[case_samples], control_df.set_index("gene_name").loc[merged_df.index, control_samples], axis=1
# ).pvalue

# # Compute ranking metric: Signed LFC weighted by -log10(p-value)
# ranking_metric = log_fc * -np.log10(ttest_pvalues)

# # Create DataFrame for ranked genes
# ranked_genes = pd.DataFrame({"gene": merged_df.index, "rank": ranking_metric})
# ranked_genes = ranked_genes.sort_values("rank", ascending=False)

# # Save ranked gene list for GSEA
# ranked_genes.to_csv("ranked_gene_list.rnk", sep="\t", index=False, header=False)

# print("Top ranked genes for GSEA:")
# print(ranked_genes.head(10))
