import pandas as pd

def standardize_patient_id(patient_id):
    """Convert TCGA ID format from TCGA.XX.XXXX to TCGA-XX-XXXX"""
    if isinstance(patient_id, str):
        return patient_id.replace('.', '-')
    return patient_id

def merge_patient_data(survival_file, split_file):
    """
    Merge patient survival data with train/validation split information.
    
    Parameters:
    survival_file (str): Path to survival data CSV
    split_file (str): Path to split information CSV
    
    Returns:
    pd.DataFrame: Merged dataset with split information
    """
    # Read the data files
    survival_data = pd.read_csv(survival_file)
    split_data = pd.read_csv(split_file)
    
    # Standardize patient IDs in both dataframes
    split_data['patient_id'] = split_data['patient_id'].apply(standardize_patient_id)
    
    # Create a dictionary from split data for faster lookup
    split_dict = dict(zip(split_data['patient_id'], split_data['split']))
    
    # Add split information to survival data
    survival_data['split'] = survival_data['patient_barcode'].map(split_dict)
    
    # Sort by patient barcode
    survival_data = survival_data.sort_values('patient_barcode')
    
    return survival_data

# Example usage:
survival_file = '/Users/stanleychen/git/Melanoma/data/patient_survival_analysis.csv'
split_file = '/Users/stanleychen/git/Melanoma/data/training_validation_split.csv'

merged_data = merge_patient_data(survival_file, split_file)

# Save the merged data
output_file = 'patient_survival_with_splits.csv'
merged_data.to_csv(output_file, index=False)

# Print sample of results
print("\nFirst few rows of merged data:")
print(merged_data.head())

# Print summary statistics
print("\nSplit distribution:")
print(merged_data['split'].value_counts())

print("\nNumber of patients without split information:")
print(merged_data['split'].isna().sum())