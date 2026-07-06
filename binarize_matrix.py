import pandas as pd
import numpy as np

def binarize_csv(input_file, output_file):
    """
    Convert all nonzero values in a CSV file to 1, while preserving zero values.
    
    Parameters:
    input_file (str): Path to the input CSV file
    output_file (str): Path to save the output CSV file
    """
    try:
        # Read the CSV file
        df = pd.read_csv(input_file)
        
        # Get the first column name (assuming it's the identifier column)
        id_column = df.columns[0]
        
        # Store the identifier column
        identifiers = df[id_column]
        
        # Process all columns except the first one (identifier column)
        for column in df.columns[1:]:
            # Convert column to numeric, coercing errors to NaN
            df[column] = pd.to_numeric(df[column], errors='coerce')
            
            # Replace all nonzero values with 1
            # Note: This will also convert NaN to 0
            df[column] = df[column].apply(lambda x: 1 if x != 0 and not pd.isna(x) else 0)
        
        # Save the processed dataframe
        df.to_csv(output_file, index=False)
        print(f"Successfully processed file. Output saved to: {output_file}")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

# Set your input and output file paths here
input_file = "/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix.csv"  # Replace with your input file path
output_file = "/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix_binary.csv"  # Replace with your desired output file path

# Run the conversion
binarize_csv(input_file, output_file)