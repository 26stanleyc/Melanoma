import pandas as pd

def filter_selected_mirnas(mirna_file, selected_file, output_file):
    """
    Filter miRNA expression data to keep only selected miRNAs.
    
    Parameters:
    -----------
    mirna_file : str
        Path to the normalized combined miRNA expression file
    selected_file : str
        Path to the file containing selected miRNA features
    output_file : str
        Path where the filtered data will be saved
    """
    # Read files
    print("Reading files...")
    mirna_df = pd.read_csv(mirna_file, index_col=0).T  # Transpose to get miRNAs as columns
    selected_df = pd.read_csv(selected_file)
    
    # Get list of selected miRNAs
    selected_mirnas = selected_df['Feature'].tolist()
    
    # Check which selected miRNAs are in the expression data
    available_mirnas = [mirna for mirna in selected_mirnas if mirna in mirna_df.columns]
    missing_mirnas = [mirna for mirna in selected_mirnas if mirna not in mirna_df.columns]
    
    if missing_mirnas:
        print("\nWARNING: The following selected miRNAs were not found in the expression data:")
        for mirna in missing_mirnas:
            print(f"- {mirna}")
    
    # Select only available miRNAs
    filtered_df = mirna_df[available_mirnas]
    
    # Save filtered dataset
    filtered_df.to_csv(output_file)
    
    print(f"\nSummary:")
    print(f"Original number of miRNAs: {len(mirna_df.columns)}")
    print(f"Number of selected miRNAs found: {len(available_mirnas)}")
    print(f"Number of selected miRNAs missing: {len(missing_mirnas)}")
    print(f"Filtered data saved to: {output_file}")
    
    # Print first few rows of filtered data for verification
    print("\nFirst few rows of filtered data:")
    print(filtered_df.head())
    
    return filtered_df

if __name__ == "__main__":
    # File paths
    mirna_file = "/Users/stanleychen/git/Melanoma/apaper/miRNA.csv"
    selected_file = "/Users/stanleychen/git/Melanoma/data/selected_features.csv"
    output_file = "/Users/stanleychen/git/Melanoma/apaper/miRNA_IBCGA.csv"
    
    print("Starting filtering process...")
    filtered_df = filter_selected_mirnas(mirna_file, selected_file, output_file)