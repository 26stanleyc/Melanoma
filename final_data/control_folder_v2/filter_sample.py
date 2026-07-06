import pandas as pd

def filter_skin_samples(phenotype_path, expression_path):
   # Read the phenotype data 
   phenotype_df = pd.read_csv(phenotype_path, sep='\t')
   
   # Get list of sample IDs that are skin tissue
   skin_samples = phenotype_df[phenotype_df['_primary_site'] == 'Skin']['Sample'].tolist()
   
   # Read expression data
   expression_df = pd.read_csv(expression_path, sep='\t')
   
   # Get all column names that have any of the skin sample IDs in them
   # We keep first column (gene IDs) and any columns containing skin sample IDs
   cols_to_keep = ['sample']  # Keep first column with gene IDs
   for col in expression_df.columns[1:]:  # Skip first column when checking
       if any(skin_id in col for skin_id in skin_samples):
           cols_to_keep.append(col)
           
   # Filter expression data to only keep relevant columns
   filtered_df = expression_df[cols_to_keep]
   
   return filtered_df

# Example usage:
filtered_data = filter_skin_samples('/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/GTEX_phenotype', '/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/gtex_RSEM_gene_fpkm')
filtered_data.to_csv('filtered_expression.txt', sep='\t', index=False)