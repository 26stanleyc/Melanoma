import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

def normalize_expression_data(data_path, output_path=None):
    """
    Normalize expression data using log2(x+1) transformation and min-max scaling
    
    Args:
        data_path (str): Path to input expression data file
        output_path (str, optional): Path to save normalized data
    
    Returns:
        pd.DataFrame: Normalized expression data
    """
    # Load data
    print(f"Loading data from {data_path}")
    df = pd.read_csv(data_path, index_col=0)
    print(f"Initial shape: {df.shape}")
    
    # 1. Apply log2(x+1) transformation
    print("Applying log2(x+1) transformation...")
    # df_log = np.log2(df + 1)
    df_log=df
    
    # 2. Apply min-max scaling to get values between 0 and 1
    print("Applying min-max normalization...")
    scaler = MinMaxScaler()
    normalized_values = scaler.fit_transform(df_log)
    
    # Convert back to DataFrame with original indices and columns
    df_normalized = pd.DataFrame(
        normalized_values,
        index=df.index,
        columns=df.columns
    )
    
    # Print value ranges to verify normalization
    print("\nValue ranges after normalization:")
    print(f"Min value: {df_normalized.values.min():.4f}")
    print(f"Max value: {df_normalized.values.max():.4f}")
    
    # Save normalized data if output path is provided
    if output_path:
        df_normalized.to_csv(output_path)
        print(f"Saved normalized data to {output_path}")
    
    return df_normalized

# Example usage
if __name__ == "__main__":
    # Replace these paths with your actual file paths
    input_path = "/Users/stanleychen/git/Melanoma/apaper/miRNA.csv"
    output_path = "/Users/stanleychen/git/Melanoma/apaper/miRNA.csv"
    
    normalized_data = normalize_expression_data(
        data_path=input_path,
        output_path=output_path
    )