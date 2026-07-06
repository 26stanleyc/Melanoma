import pandas as pd
from scipy import stats

# File paths (replace with actual file paths if needed)
gene_expression_file = "/Users/stanleychen/git/Melanoma/final_data/cleaned_RNA-seq_data.csv"
label_file = "/Users/stanleychen/git/Melanoma/final_data/survclass_filtered_patient_data_all_v2.csv"

# Load the gene expression data
gene_data = pd.read_csv(gene_expression_file, index_col=0)  # Set first column (genes) as index

# Transpose so that 'patient_barcode' becomes a column
gene_data = gene_data.T
gene_data.index.name = "patient_barcode"  # Rename index to match the label file

# Load the survival label data
label_data = pd.read_csv(label_file, usecols=['patient_barcode', 'survival_class'])

# Merge gene expression data with survival labels based on 'patient_barcode'
merged_data = gene_data.merge(label_data, on='patient_barcode')

# Extract gene expression values and survival class
genes = merged_data.columns.difference(['patient_barcode', 'survival_class'])  # Select only gene columns
labels = merged_data['survival_class']

# Perform ANOVA for each gene
anova_results = {}
for gene in genes:
    groups = [merged_data[gene][labels == label].dropna() for label in labels.unique()]
    if all(len(group) > 1 for group in groups):  # Ensure valid ANOVA input
        f_val, p_val = stats.f_oneway(*groups)
        anova_results[gene] = {'F-Value': f_val, 'P-Value': p_val}

# Convert results to a DataFrame
anova_results_df = pd.DataFrame.from_dict(anova_results, orient='index')

# Sort by F-Value to identify the most significant genes
anova_results_df = anova_results_df.sort_values(by='F-Value', ascending=False)

# Select the top 500 genes
top_500_genes = anova_results_df.head(500)

# Save or use the results
top_500_genes.to_csv("top_500_anova_genes.csv")

# Store in a variable for further analysis
imported_anova_results = top_500_genes
