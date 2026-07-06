import pandas as pd

def count_unique_genes(data_file):
    # Read the data file
    df = pd.read_csv(data_file)
    
    # Count unique genes
    unique_genes = df['gene'].nunique()
    
    # Get the list of unique genes
    gene_list = sorted(df['gene'].unique())
    
    print(f"Number of unique genes: {unique_genes}")
    # print("\nList of unique genes:")
    # for gene in gene_list:
    #     print(gene)

# Example usage (assuming data is saved in variants.csv)
if __name__ == "__main__":
    count_unique_genes("/Users/stanleychen/git/Melanoma/final_data/mutations_with_alphamissense.csv")