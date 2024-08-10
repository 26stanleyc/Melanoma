import pandas as pd

def normalize_dataframe(df):
    """
    Normalize each row of the DataFrame by subtracting the row average.
    
    Args:
    df (pandas.DataFrame): Input DataFrame
    
    Returns:
    pandas.DataFrame: Normalized DataFrame
    """
    return df.sub(df.mean(axis=1), axis=0)

# Path to your existing CSV file
input_file = '/Users/stanleychen/git/Melanoma/data/combined_rna_seq_data.csv'

# Read the CSV file
# Assuming the first column is the gene names and should be used as index
df = pd.read_csv(input_file, index_col=0)

# Display the first few rows of the original DataFrame
print("Original DataFrame:")
print(df.head())

# Normalize the DataFrame
normalized_df = normalize_dataframe(df)

# Display the first few rows of the normalized DataFrame
print("\nNormalized DataFrame:")
print(normalized_df.head())

# Save the normalized DataFrame to a new CSV file
output_file = 'normalized_rna_seq_data.csv'
normalized_df.to_csv(output_file)

print(f"\nNormalized data saved to {output_file}")