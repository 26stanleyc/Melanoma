import pandas as pd

# Load the primary data
primary_df = pd.read_csv('/Users/stanleychen/git/Melanoma/final_data/RNA_cancer_related_genes.csv', index_col=0)

# Load the mapping file, skipping the first 4 rows
mapping_df = pd.read_csv('/Users/stanleychen/git/Melanoma/clinical_data/data/RNA-seq/0aaec814-e30d-489b-928d-62776b8028cd/dd28c353-13f2-4e6a-a2d7-df94352e043f.rna_seq.augmented_star_gene_counts.tsv', 
                        sep='\t',
                        skiprows=4)

# Clean up the gene IDs by removing version numbers (e.g., ".15" from "ENSG00000000003.15")
mapping_df.iloc[:, 0] = mapping_df.iloc[:, 0].str.split('.').str[0]

# Create a dictionary for mapping gene names to ENSG IDs
gene_to_ensg = dict(zip(mapping_df['Unnamed: 1'], mapping_df.iloc[:, 0]))

# Map the index (gene names) of the primary DataFrame to ENSG IDs
primary_df.index = primary_df.index.map(gene_to_ensg)

# Drop rows (genes) that couldn't be mapped
primary_df.dropna(inplace=True)

# Transpose the DataFrame so that genes become columns
transposed_df = primary_df.transpose()

# Reset index to turn the patient identifiers into a column
transposed_df.reset_index(inplace=True)

# Rename the index column to 'identifier'
transposed_df.rename(columns={'index': 'identifier'}, inplace=True)

# Move the 'identifier' column to the end
identifier_col = transposed_df.pop('identifier')
transposed_df['identifier'] = identifier_col

# Save the transformed DataFrame to a new CSV file
transposed_df.to_csv('transformed_data.csv', index=False)