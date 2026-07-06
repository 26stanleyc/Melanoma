import pandas as pd
import numpy as np

def convert_nans_to_zero(input_file, output_file):
    """
    Read CSV file, convert NaN values to 0, and save to new file
    """
    # Read the CSV file
    print(f"Reading file: {input_file}")
    df = pd.read_csv(input_file)
    
    # Print initial info
    print("\nInitial data info:")
    print(f"Shape: {df.shape}")
    print(f"Number of NaN values: {df.isna().sum().sum()}")
    
    # Convert NaN to 0
    df = df.fillna(0)
    
    # Verify conversion
    print("\nAfter conversion:")
    print(f"Number of remaining NaN values: {df.isna().sum().sum()}")
    
    # Save to new file
    df.to_csv(output_file, index=False)
    print(f"\nSaved cleaned data to: {output_file}")
    
    return df

if __name__ == "__main__":
    input_file = "/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix.csv"  # Replace with your input file name
    output_file = "/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix_with_0s.csv"    # Replace with desired output file name
    
    cleaned_df = convert_nans_to_zero(input_file, output_file)