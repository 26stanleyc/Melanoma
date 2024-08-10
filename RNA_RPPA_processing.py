import pandas as pd
import numpy as np
from pathlib import Path
import os

def fpkm_to_tpm(fpkm):
    return fpkm / (np.sum(fpkm) / 1e6)

def process_gene_expression(input_file):
    # Read the gene expression data
    df = pd.read_csv(input_file, sep='\t', comment='#')
    
    # Select relevant columns
    df = df[['gene_name', 'unstranded']]
    df.columns = ['gene', 'fpkm']
    
    # Filter out non-protein coding genes
    df = df[df['gene'] != 'protein_coding']
    
    # Convert FPKM to TPM
    df['tpm'] = fpkm_to_tpm(df['fpkm'])
    
    # Apply log2(TPM + 1)
    df['log2_tpm_plus_1'] = np.log2(df['tpm'] + 1)
    
    # Keep only the relevant columns
    df = df[['gene', 'log2_tpm_plus_1']]
    
    return df

def process_protein_expression(input_file):
    # Read the protein expression data
    df = pd.read_csv(input_file, sep='\t')
    
    # Select relevant columns
    df = df[['peptide_target', 'protein_expression']]
    
    # Set peptide_target as index
    df.set_index('peptide_target', inplace=True)
    
    # Remove rows with NA values
    df.dropna(inplace=True)
    
    return df

def process_file(input_file):
    input_path = Path(input_file)
    
    if 'gene' in input_file.lower() or 'rna' in input_file.lower():
        df = process_gene_expression(input_file)
    elif 'protein' in input_file.lower() or 'rppa' in input_file.lower():
        df = process_protein_expression(input_file)
    else:
        print(f"Unrecognized file type: {input_file}")
        return
    
    # Create output filename
    if 'gene' in input_file.lower() or 'rna' in input_file.lower():
        output_file = Path('/Users/stanleychen/git/Melanoma/data/RNA-seq_processed/a').with_name(f"{input_path.stem}_processed{input_path.suffix}")
    else:
        output_file = Path('/Users/stanleychen/git/Melanoma/data/RPPA_processed/a').with_name(f"{input_path.stem}_processed{input_path.suffix}")
    
    # Save processed data
    df.to_csv(output_file, index=False)  # Save without index to keep the format consistent
    print(f"Processed data saved to {output_file}")

# Example usage
root_dir='data/RNA-seq'

for root, dirs, files in os.walk(root_dir):
    # Skip the root directory itself
    if root == root_dir:
        continue
    # Skip 'logs' directories
    if 'logs' in root:
        continue
    
    for file in files:
        if file.endswith('.tsv'):  # Only consider TSV files
            file_path = os.path.join(root, file)
            # print(file_path + '\t')
            process_file(file_path)

# Uncomment and modify the following block to process RPPA data if needed
# root_dir='RPPA'
# for root, dirs, files in os.walk(root_dir):
#     # Skip the root directory itself
#     if root == root_dir:
#         continue
#     # Skip 'logs' directories
#     if 'logs' in root:
#         continue
    
#     for file in files:
#         if file.endswith('.tsv'):  # Only consider TSV files
#             file_path = os.path.join(root, file)
#             # print(file_path + '\t')
#             process_file(file_path)
