import pandas as pd
import numpy as np

def aggregate_patient_data(group):
    result = {}
    for column in group.columns:
        if column in ['bcr_patient_uuid', 'bcr_patient_barcode']:
            result[column] = group[column].iloc[0]
        elif pd.api.types.is_numeric_dtype(group[column]):
            result[column] = group[column].mean()
        elif pd.api.types.is_datetime64_any_dtype(group[column]):
            result[column] = group[column].max()
        else:
            # For categorical data, concatenate unique values
            unique_values = group[column].dropna().unique()
            result[column] = ', '.join(unique_values) if len(unique_values) > 0 else np.nan
    return pd.Series(result)

def process_melanoma_data(data_files):
    dataframes = {}
    
    for file in data_files:
        df = pd.read_csv(file, sep='\t')
        file_prefix = file.split('_')[1]
        
        if 'bcr_patient_uuid' not in df.columns or 'bcr_patient_barcode' not in df.columns:
            print(f"Error: {file} is missing required identifier columns.")
            continue
        
        df.columns = [f"{file_prefix}_{col}" if col not in ['bcr_patient_uuid', 'bcr_patient_barcode'] else col for col in df.columns]
        dataframes[file_prefix] = df
        print(f"Loaded {file}: {df.shape[0]} rows, {df.shape[1]} columns")
    
    # Merge all dataframes
    merged_data = None
    for df in dataframes.values():
        if merged_data is None:
            merged_data = df
        else:
            merged_data = pd.merge(merged_data, df, on=['bcr_patient_uuid', 'bcr_patient_barcode'], how='outer')
    
    print(f"Merged data shape before deduplication: {merged_data.shape}")
    
    # Deduplicate and aggregate data
    deduplicated_data = merged_data.groupby(['bcr_patient_uuid', 'bcr_patient_barcode']).apply(aggregate_patient_data).reset_index(drop=True)
    
    print(f"Deduplicated data shape: {deduplicated_data.shape}")
    
    return deduplicated_data

# Usage
data_files = [
    'clinical_patient_processed.tsv', 
    'clinical_drug_processed.tsv', 
    'clinical_radiation_processed.tsv',
    'clinical_followup_processed.tsv',
    'clinical_omf_processed.tsv',
    'clinical_tme_processed.tsv'
]

combined_data = process_melanoma_data(data_files)

print("\nFinal Combined Data Shape:", combined_data.shape)
print("Sample columns in final data:", list(combined_data.columns)[:20])

# Save the combined data to a file
combined_data.to_csv('combined_clinical_data.tsv', sep='\t', index=False)

# Print a sample of the deduplicated data
print("\nSample of deduplicated data:")
print(combined_data.head())

# Check for any remaining duplicates
duplicate_check = combined_data.duplicated(subset=['bcr_patient_uuid', 'bcr_patient_barcode'])
if duplicate_check.any():
    print(f"\nWarning: {duplicate_check.sum()} duplicate entries found after deduplication.")
else:
    print("\nNo duplicate entries found after deduplication.")