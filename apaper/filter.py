import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer

def preprocess_omics_data(data_path, output_path=None, sample_threshold=0.5, zero_threshold=0.2):
    """
    Preprocess omics data according to specified rules:
    1. Remove probes with > 50% missing values
    2. Impute remaining missing values using KNN
    3. Remove genes with > 20% zero values
    
    Args:
        data_path (str): Path to input data file
        output_path (str, optional): Path to save processed data
        sample_threshold (float): Threshold for missing values (default 0.5)
        zero_threshold (float): Threshold for zero values (default 0.2)
    
    Returns:
        pd.DataFrame: Processed data
    """
    # Load data
    print(f"Loading data from {data_path}")
    df = pd.read_csv(data_path, index_col=0)
    print(f"Initial shape: {df.shape}")
    
    # 1. Remove probes with > 50% missing values
    missing_counts = df.isna().sum(axis=1)
    missing_mask = missing_counts <= (df.shape[1] * sample_threshold)
    df = df[missing_mask]
    print(f"Removed {sum(~missing_mask)} probes with >{sample_threshold*100}% missing values")
    print(f"Shape after removing high-missing probes: {df.shape}")
    
    # 2. KNN imputation
    print("Performing KNN imputation...")
    imputer = KNNImputer(n_neighbors=5)
    df_imputed = pd.DataFrame(
        imputer.fit_transform(df),
        index=df.index,
        columns=df.columns
    )
    
    # 3. Remove genes (rows) with > 20% zeros
    zero_counts_per_gene = (df_imputed == 0).sum(axis=1)  # Count zeros per row/gene
    zero_percentage_per_gene = zero_counts_per_gene / df_imputed.shape[1]
    genes_to_keep = zero_percentage_per_gene <= zero_threshold
    
    df_final = df_imputed.loc[genes_to_keep, :]
    
    print(f"Removed {sum(~genes_to_keep)} genes with >{zero_threshold*100}% zeros")
    print(f"Final shape: {df_final.shape}")
    
    # Save processed data if output path is provided
    if output_path:
        df_final.to_csv(output_path)
        print(f"Saved processed data to {output_path}")
    
    return df_final

# Example usage
if __name__ == "__main__":
    # Replace these paths with your actual file paths
    input_path = "/Users/stanleychen/git/Melanoma/data/combined_miRNA.csv"
    output_path = "/Users/stanleychen/git/Melanoma/apaper/miRNA.csv"
    
    # Process the data
    processed_data = preprocess_omics_data(
        data_path=input_path,
        output_path=output_path,
        sample_threshold=0.5,  # Remove features with >50% missing values
        zero_threshold=0.2     # Remove genes with >20% zeros
    )