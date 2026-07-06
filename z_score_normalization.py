import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

def read_data_file(file_path: str) -> pd.DataFrame:
    """
    Read data from CSV or TSV file
    
    Args:
        file_path: Path to the input file
        
    Returns:
        pandas DataFrame with gene expression data
    """
    # Get file extension
    file_ext = Path(file_path).suffix.lower()
    
    try:
        if file_ext == '.csv':
            df = pd.read_csv(file_path)
        elif file_ext == '.tsv':
            df = pd.read_csv(file_path, sep='\t')
        else:
            raise ValueError("File must be either .csv or .tsv format")
            
        # Ensure first column is named 'gene'
        df.columns.values[0] = 'gene'
        return df
        
    except Exception as e:
        print(f"Error reading file {file_path}: {str(e)}")
        raise

def zscore_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform Z-score normalization on the numeric columns of the dataframe
    
    Args:
        df: Input DataFrame with 'gene' column and numeric data columns
        
    Returns:
        Normalized DataFrame
    """
    # Separate gene names and numeric data
    genes = df['gene']
    numeric_data = df.select_dtypes(include=[np.number])
    
    # Perform Z-score normalization
    normalized_data = (numeric_data - numeric_data.mean()) / numeric_data.std()
    
    # Add gene names back
    normalized_data.insert(0, 'gene', genes)
    
    return normalized_data

def get_summary_stats(original_df: pd.DataFrame, normalized_df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """
    Calculate summary statistics for original and normalized data
    
    Args:
        original_df: Original DataFrame
        normalized_df: Normalized DataFrame
        dataset_name: Name of the dataset
        
    Returns:
        DataFrame with summary statistics
    """
    numeric_orig = original_df.select_dtypes(include=[np.number])
    numeric_norm = normalized_df.select_dtypes(include=[np.number])
    
    stats = {
        'Dataset': dataset_name,
        'Original_Mean': numeric_orig.values.mean(),
        'Original_Std': numeric_orig.values.std(),
        'Normalized_Mean': numeric_norm.values.mean(),
        'Normalized_Std': numeric_norm.values.std()
    }
    
    return pd.DataFrame([stats])

def plot_distributions(original_df: pd.DataFrame, normalized_df: pd.DataFrame, 
                      output_dir: str, dataset_name: str):
    """
    Create distribution plots before and after normalization
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot original data distribution
    orig_data = original_df.select_dtypes(include=[np.number]).values.flatten()
    ax1.hist(orig_data, bins=50)
    ax1.set_title(f'Original Data Distribution\n{dataset_name}')
    ax1.set_xlabel('Expression Values')
    ax1.set_ylabel('Frequency')
    
    # Plot normalized data distribution
    norm_data = normalized_df.select_dtypes(include=[np.number]).values.flatten()
    ax2.hist(norm_data, bins=50)
    ax2.set_title(f'Normalized Data Distribution\n{dataset_name}')
    ax2.set_xlabel('Z-score')
    ax2.set_ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{dataset_name}_distribution.png')
    plt.close()

def normalize_expression_data(file_path1: str, file_path2: str, output_dir: str = 'output'):
    """
    Main function to normalize gene expression data
    
    Args:
        file_path1: Path to first input file
        file_path2: Path to second input file
        output_dir: Directory for output files
    """
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Read input files
    print("Reading input files...")
    df1 = read_data_file(file_path1)
    df2 = read_data_file(file_path2)
    
    # Normalize data
    print("Performing Z-score normalization...")
    norm_df1 = zscore_normalize(df1)
    norm_df2 = zscore_normalize(df2)
    
    # Get summary statistics
    stats1 = get_summary_stats(df1, norm_df1, "Dataset1")
    stats2 = get_summary_stats(df2, norm_df2, "Dataset2")
    summary_stats = pd.concat([stats1, stats2])
    
    # Create plots
    print("Creating distribution plots...")
    plot_distributions(df1, norm_df1, output_dir, "Dataset1")
    plot_distributions(df2, norm_df2, output_dir, "Dataset2")
    
    # Save results
    print("Saving normalized data and statistics...")
    norm_df1.to_csv(f'{output_dir}/normalized_dataset1.csv', index=False)
    norm_df2.to_csv(f'{output_dir}/normalized_dataset2.csv', index=False)
    summary_stats.to_csv(f'{output_dir}/normalization_statistics.csv', index=False)
    
    print(f"Results saved to {output_dir}")
    print("\nSummary Statistics:")
    print(summary_stats)

# Example usage:
if __name__ == "__main__":
    normalize_expression_data(
        "/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/filtered_rna_seq-R.csv",
        "/Users/stanleychen/git/Melanoma/final_data/control_folder/filtered_control-R.tsv",
        "/Users/stanleychen/git/Melanoma/final_data/data_distribution"
    )