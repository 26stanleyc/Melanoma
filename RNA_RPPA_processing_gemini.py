import pandas as pd
import numpy as np
from pathlib import Path
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_gene_expression(input_file):
    """Processes gene expression data (TPM) and applies log2 transformation.

    Args:
        input_file (str): Path to the gene expression file (TSV).

    Returns:
        pandas.DataFrame: A DataFrame with log2(TPM + 1) values, or None if an error occurs.
    """
    try:
        # Read the gene expression data, handling potential comment lines
        df = pd.read_csv(input_file, sep='\t', comment='#')

        # Select relevant columns. Adjust column names if different.
        if 'tpm_unstranded' in df.columns:
            tpm_col = 'tpm_unstranded'
        elif 'TPM' in df.columns:
            tpm_col = 'TPM'
        else:
            raise KeyError("TPM column not found. Check input file format")

        if 'gene_name' in df.columns:
            gene_name_col = 'gene_name'
        else:
            raise KeyError("gene_name column not found. Check input file format")
        df = df[[gene_name_col, tpm_col]]
        df.columns = ['gene', 'tpm']

        # Convert TPM column to numeric, handling potential errors
        df['tpm'] = pd.to_numeric(df['tpm'], errors='coerce')
        df = df.dropna(subset=['tpm']) #remove rows with NaN values after conversion
        df = df[df['tpm'] > 0] #added this

        # Apply log2(TPM + 1)
        df['log2_tpm_plus_1'] = np.log2(df['tpm'] + 1)

        # Set 'gene' as index AFTER processing
        df.set_index('gene', inplace=True)

        return df

    except FileNotFoundError:
        logging.error(f"Error: Input file '{input_file}' not found.")
        return None
    except KeyError as e:
        logging.error(f"Error: Column '{e}' not found in the input file. Check column names.")
        return None
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        return None

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
        output_dir = '/Users/stanleychen/git/Melanoma/data/RNA-seq_processed_gemini/a'
    elif 'protein' in input_file.lower() or 'rppa' in input_file.lower():
        # df = process_protein_expression(input_file)
        # output_dir = '/Users/stanleychen/git/Melanoma/data/RPPA_processed/a'
        print('protein')
    else:
        logging.warning(f"Unrecognized file type: {input_file}")
        return
    
    if df is None: #added this
        return

    os.makedirs(output_dir, exist_ok=True) #added this
    # Create output filename
    output_file = Path(output_dir).with_name(f"{input_path.stem}_processed{input_path.suffix}")
    
    # Save processed data
    df.to_csv(output_file, index=True)  # Save with index
    logging.info(f"Processed data saved to {output_file}")

# Example usage
root_dir='RNA-seq'

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
root_dir='RPPA'
for root, dirs, files in os.walk(root_dir):
    # Skip the root directory itself
    if root == root_dir:
        continue
    # Skip 'logs' in root:
    if 'logs' in root:
        continue
    
    for file in files:
        if file.endswith('.tsv'):  # Only consider TSV files
            file_path = os.path.join(root, file)
            # print(file_path + '\t')
            process_file(file_path)