import pandas as pd
import numpy as np
from pathlib import Path
import random

def sample_columns(input_file, n_samples=500):
    # Read the existing processed file
    df = pd.read_csv(input_file, sep='\t')
    
    # Keep gene column and sample n random columns
    sample_cols = random.sample(list(df.columns[1:]), min(n_samples, len(df.columns)-1))
    selected_cols = ['gene'] + sample_cols
    
    # Create new dataframe with selected columns
    result_df = df[selected_cols]
    
    # Save to new file
    output_file = str(input_file).replace('.tsv', f'_{n_samples}samples.tsv')
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved file with {len(result_df)} genes and {len(sample_cols)} samples")

# Example usage
input_file = "/Users/stanleychen/git/Melanoma/final_data/control_folder/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm_processed.tsv"
sample_columns(input_file, 500)