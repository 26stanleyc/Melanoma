import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

def process_gtex_expression(input_file):
    print("Reading GCT file...")
    
    with open(input_file, 'r') as f:
        next(f)
        dims = next(f).strip().split('\t')
        n_genes, n_samples = int(dims[0]), int(dims[1])
        print(f"Found {n_genes} genes and {n_samples} samples")
    
    df = pd.read_csv(input_file, sep='\t', skiprows=2)
    gene_cols = df.columns[:2].tolist()
    sample_cols = df.columns[2:].tolist()
    
    print(f"Processing {len(sample_cols)} samples...")
    
    # Pre-allocate the log2 transformed data matrix
    log2_data = np.log2(df[sample_cols].values + 1)
    
    # Create result dataframe all at once
    result_df = pd.DataFrame(
        log2_data,
        columns=sample_cols,
        index=df[gene_cols[1]]
    ).reset_index()
    result_df.columns = ['gene'] + sample_cols
    
    return result_df

def process_file(input_file):
    input_path = Path(input_file)
    print(f"\nProcessing file: {input_path.name}")
    print("=" * 50)
    
    df = process_gtex_expression(input_file)
    
    output_file = input_path.parent / f"{input_path.stem}_processed.tsv"
    
    print(f"\nSaving processed data to {output_file}")
    with tqdm(total=1, desc="Saving file") as pbar:
        df.to_csv(output_file, sep='\t', index=False)
        pbar.update(1)
    
    print("\nProcessing Summary:")
    print(f"Total genes processed: {len(df)}")
    print(f"Total samples: {len(df.columns) - 1}")
    
    print("\nFirst few rows of processed data:")
    print(df.head())

def main():
    input_file = "/Users/stanleychen/git/Melanoma/final_data/control_folder/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct"
    
    try:
        process_file(input_file)
        print("\nProcessing completed successfully!")
    except Exception as e:
        print(f"\nError during processing: {str(e)}")
        raise e

if __name__ == "__main__":
    main()