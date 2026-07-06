import pandas as pd
import numpy as np

def filter_zero_rows(input_file, output_file, zero_threshold=0.3):
    """
    Filter out rows from a CSV file where more than a specified threshold of columns contain zero values.
    
    Parameters:
    -----------
    input_file : str
        Path to the input CSV file
    output_file : str
        Path where the filtered CSV will be saved
    zero_threshold : float, optional
        Maximum allowed proportion of zero values in a row (default: 0.3)
        
    Returns:
    --------
    tuple
        (number of rows before filtering, number of rows after filtering)
    """
    try:
        # Read the CSV file
        df = pd.read_csv(input_file, index_col=0)
        
        # Store original number of rows
        original_rows = len(df)
        
        # Calculate the proportion of zeros in each row
        zero_proportions = (df == 0).sum(axis=1) / df.shape[1]
        
        # Keep only rows where the proportion of zeros is <= threshold
        df_filtered = df[zero_proportions <= zero_threshold]
        
        # Save the filtered dataset
        df_filtered.to_csv(output_file)
        
        return original_rows, len(df_filtered)
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        return None, None

if __name__ == "__main__":
    # Example usage
    input_file = "/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA.csv"  # Change this to your input file name
    output_file = "/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA.csv"  # Change this to your desired output file name
    
    before, after = filter_zero_rows(input_file, output_file)
    
    if before is not None and after is not None:
        print(f"Original number of rows: {before}")
        print(f"Rows after filtering: {after}")
        print(f"Removed {before - after} rows with more than 30% zero values")