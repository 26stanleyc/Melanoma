import pandas as pd

def get_unique_data_size(file_path):
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Print column names to see what we're working with
    print("Available columns:", df.columns.tolist())
    
    # Get the first column name (assuming gene names are in the first column)
    gene_column = df.columns[0]
    
    # Drop duplicates based on the first column
    df_unique = df.drop_duplicates(subset=[gene_column])
    
    return len(df_unique)

data = "/Users/stanleychen/git/Melanoma/final_data/control_folder/control_sample_processed.tsv"
print("Number of unique entries:", get_unique_data_size(data))