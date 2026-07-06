import pandas as pd
import numpy as np

def clean_and_filter_matrix(input_path, output_path, zero_threshold=0.3):
    """
    Clean matrix by converting NaN to 0, then remove rows with too many zeros.
    
    Args:
        input_path (str): Path to input CSV file
        output_path (str): Path to save filtered CSV
        zero_threshold (float): Maximum allowed proportion of zeros (default 0.3)
    """
    # Read the CSV file
    print("Reading input file...")
    df = pd.read_csv(input_path, index_col=0)
    
    initial_rows = df.shape[0]
    print(f"Initial number of rows: {initial_rows}")
    
    # Count NaN values before conversion
    nan_counts = df.isna().sum().sum()
    print(f"\nFound {nan_counts} NaN values in the matrix")
    
    # Convert NaN to 0
    df = df.fillna(0)
    print("Converted all NaN values to 0")
    
    # Calculate proportion of zeros in each row
    # Consider both exact zeros and very small values close to zero
    zero_proportions = ((df == 0) | (df.abs() < 1e-10)).mean(axis=1)
    
    # Keep rows with fewer zeros than threshold
    df_filtered = df[zero_proportions < zero_threshold]
    
    remaining_rows = df_filtered.shape[0]
    removed_rows = initial_rows - remaining_rows
    
    print(f"\nFiltering Results:")
    print(f"Rows removed: {removed_rows}")
    print(f"Rows remaining: {remaining_rows}")
    print(f"Proportion removed: {(removed_rows/initial_rows)*100:.2f}%")
    
    # Save filtered matrix
    print(f"\nSaving filtered matrix to {output_path}")
    df_filtered.to_csv(output_path)
    
    # Print some example rows that were removed
    print("\nExample rows removed (showing first 5):")
    removed_rows_df = df[zero_proportions >= zero_threshold]
    for idx, row in removed_rows_df.head().iterrows():
        zero_count = ((row == 0) | (row.abs() < 1e-10)).sum()
        print(f"Gene {idx}: {zero_count} zeros out of {len(row)} values ({zero_count/len(row)*100:.1f}%)")
        
    # Print summary of remaining data
    print("\nSummary of filtered matrix:")
    print(f"Shape: {df_filtered.shape}")
    print("\nSample of non-zero values distribution in first 5 rows:")
    print(df_filtered.head().describe())

# File paths
INPUT_PATH = "/Users/stanleychen/git/Melanoma/final_data/combined_RNA-seq_data.csv"  # Replace with your input file path
OUTPUT_PATH = "/Users/stanleychen/git/Melanoma/final_data/cleaned_RNA-seq_data.csv"   # Replace with desired output path

if __name__ == "__main__":
    clean_and_filter_matrix(INPUT_PATH, OUTPUT_PATH, zero_threshold=0.3)