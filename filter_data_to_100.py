import pandas as pd

def filter_expression_data():
    # Read the top 100 genes file
    top_genes_df = pd.read_csv('/Users/stanleychen/git/Melanoma/data/top_100_differential_expression_after_clean.csv')
    
    # Get the list of genes we want to keep
    genes_to_keep = set(top_genes_df['gene'].tolist())
    
    # Read the original expression data
    expression_data = pd.read_csv('/Users/stanleychen/git/Melanoma/data/cleaned_RNA_expression_data.csv')
    
    # Set the first column (gene names) as index
    expression_data.set_index(expression_data.columns[0], inplace=True)
    
    # Filter to keep only the genes in our top 100 list
    filtered_data = expression_data[expression_data.index.isin(genes_to_keep)]
    
    # Save the filtered data
    filtered_data.to_csv('filtered_expression_top_100_dx_after_clean.csv')
    
    print(f"Original number of genes: {len(expression_data)}")
    print(f"Number of genes after filtering: {len(filtered_data)}")

if __name__ == "__main__":
    filter_expression_data()