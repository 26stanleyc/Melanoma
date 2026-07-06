import pandas as pd
import numpy as np

def convert_to_log2_tpm(input_file):
    # Read data
    df = pd.read_csv(input_file, sep='\t')
    
    # Store gene IDs (first column)
    genes = df.iloc[:, 0]
    
    # Get numeric values
    values = df.iloc[:, 1:].values
    
    # Convert from log2(FPKM + 0.001) back to FPKM
    fpkm_values = (2**values) - 0.001
    
    # Convert negative values (due to floating point) to 0
    fpkm_values = np.maximum(fpkm_values, 0)
    
    # Convert to TPM
    # Sum each column and calculate scaling factor
    scaling_factors = fpkm_values.sum(axis=0) / 1e6
    
    # Convert to TPM
    tpm_values = fpkm_values / scaling_factors[None, :]
    
    # Convert to log2(TPM + 1)
    log2_tpm_values = np.log2(tpm_values + 1)
    
    # Create new dataframe
    log2_tpm_df = pd.DataFrame(log2_tpm_values, columns=df.columns[1:])
    log2_tpm_df.insert(0, df.columns[0], genes)
    
    return log2_tpm_df

# Example usage:

log2_tpm_df = convert_to_log2_tpm('/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/filtered_control-R.tsv')
log2_tpm_df.to_csv('filtered_expression_log2_tpm-R', sep='\t', index=False)