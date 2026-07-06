import pandas as pd
import numpy as np

def create_patient_gene_matrix(input_file, output_file):
    """
    Transform variant data into a matrix with:
    - Rows: patient IDs
    - Columns: gene names
    - Values: AlphaMissense scores
    """
    # Read the filtered variants file
    print("Reading input file...")
    df = pd.read_csv(input_file)
    
    # Print first few rows and columns for debugging
    print("\nFirst few rows of input data:")
    print(df.head())
    print("\nColumns in input data:")
    print(df.columns.tolist())
    
    # Create pivot table
    # If there are multiple variants per patient-gene combination,
    # use the mean score (you can change this to min, max, or other aggregation)
    print("\nCreating patient-gene matrix...")
    matrix = df.pivot_table(
        index='patient_barcode',
        columns='gene',
        values='AlphaMissense_Score',  # This should be the AlphaMissense score column
        aggfunc='mean',
        fill_value=np.nan
    )
    
    # Print matrix info
    print("\nMatrix shape:", matrix.shape)
    print("Number of patients:", len(matrix.index))
    print("Number of genes:", len(matrix.columns))
    
    # Save to file
    matrix.to_csv(output_file)
    print(f"\nSaved matrix to {output_file}")
    
    return matrix

def main():
    # File paths
    input_file = '/Users/stanleychen/git/Melanoma/final_data/filtered_mutation_data.csv'  # The file from the previous script
    output_file = 'patient_maf_matrix.csv'
    
    # Create matrix
    matrix = create_patient_gene_matrix(input_file, output_file)
    
    # Print sample of the matrix
    print("\nSample of the matrix (first 5 patients, first 5 genes):")
    print(matrix.iloc[:5, :5])

if __name__ == "__main__":
    main()