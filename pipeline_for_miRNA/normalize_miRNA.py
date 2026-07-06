import pandas as pd
import numpy as np

def log2_transform_rpm(input_file_path, output_file_path):
    """
    Apply log2 transformation to RPM-normalized miRNA expression data from an input CSV file
    and save the transformed data to an output CSV file.

    Parameters:
    input_file_path (str): Path to the input CSV file containing RPM values.
    output_file_path (str): Path to the output CSV file to save log2(RPM + 1) values.
    """
    # Read the input CSV file into a DataFrame
    rpm_df = pd.read_csv(input_file_path, index_col=0)

    # Apply log2 transformation with a pseudocount of 1
    log2_rpm_df = np.log2(rpm_df + 1)

    # Save the transformed DataFrame to the output CSV file
    log2_rpm_df.to_csv(output_file_path)

input_path = '/Users/stanleychen/git/Melanoma/data/combined_miRNA.csv'
output_path = '/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA.csv'

log2_transform_rpm(input_path, output_path)
