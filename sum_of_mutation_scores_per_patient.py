import pandas as pd
import numpy as np

def calculate_and_save_patient_sums(input_file, output_file):
    """
    Calculate the sum of values for each patient (row) in the matrix and save to file
    """
    # Read the CSV file
    df = pd.read_csv(input_file)
    
    # Convert empty strings and non-numeric values to NaN
    df = df.replace('', np.nan)
    
    # Get patient IDs from first column
    patient_ids = df.iloc[:, 0]
    
    # Convert remaining columns to numeric, coercing errors to NaN
    numeric_data = df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
    
    # Calculate row sums
    row_sums = numeric_data.sum(axis=1)
    
    # Create a new dataframe with just patient IDs and their sums
    results = pd.DataFrame({
        'patient_barcode': patient_ids,
        'total_score': row_sums
    })
    
    # Sort by total score descending
    results = results.sort_values('total_score', ascending=False)
    
    # Save to file without index
    results.to_csv(output_file, index=False, header=False)
    print(f"Results saved to {output_file}")
    
    return results

def main():
    input_file = '/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix.csv'  # Your input matrix file
    output_file = 'maf_sums.csv'  # Output file for results
    
    # Calculate and save results
    results = calculate_and_save_patient_sums(input_file, output_file)

if __name__ == "__main__":
    main()