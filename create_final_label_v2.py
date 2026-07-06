import pandas as pd

def convert_stage_to_grade(stage):
    """
    Convert tumor stage to numeric grade, handling edge cases
    """
    if pd.isna(stage):
        return None
        
    # Convert to string, remove whitespace, and convert to uppercase
    stage = str(stage).strip().upper()
    
    # Handle special cases
    if stage == '[NOT AVAILABLE]' or stage == 'NOT REPORTED' or stage == 'NA':
        return None
    
    if 'I/II NOS' in stage:
        return 1  # Convert I/II NOS to grade 1 (equivalent to stage II)
    
    if stage == 'STAGE 0' or stage == '0':
        return 0  # Convert stage 0 to grade 0
        
    # Check for main stages
    if 'IV' in stage:
        return 3
    elif 'III' in stage:
        return 2
    elif 'II' in stage or 'I/II' in stage:  # Including I/II cases
        return 1
    elif 'I' in stage:
        return 0
    else:
        return None

def process_patient_data(input_file, output_file):
    """
    Process patient data and create grade labels
    """
    # Read the input file
    df = pd.read_csv(input_file)
    
    # Create result dataframe with required columns
    result_df = pd.DataFrame()
    result_df['patient_barcode'] = df['patient_barcode']
    
    # Convert survival time to float and round to 1 decimal
    result_df['survival_time_days'] = pd.to_numeric(df['survival_time_days'], errors='coerce').round(1)
    
    # Add other columns
    result_df['survival_time_source'] = df['survival_time_source']
    result_df['tumor_stage'] = df['tumor_stage']  # Keep original stage for reference
    result_df['vital_status'] = df['vital_status']
    
    # Convert stage to grade
    result_df['grade'] = df['tumor_stage'].apply(convert_stage_to_grade)
    
    # Remove rows with NaN values
    initial_rows = len(result_df)
    result_df = result_df.dropna()
    removed_rows = initial_rows - len(result_df)
    
    # Sort by survival time
    result_df = result_df.sort_values('survival_time_days', ascending=False)
    
    # Save to file
    result_df.to_csv(output_file, index=False)
    
    # Print summary statistics
    print(f"\nProcessing complete:")
    print(f"Input rows: {initial_rows}")
    print(f"Rows removed: {removed_rows}")
    print(f"Final rows: {len(result_df)}")
    
    print("\nGrade distribution:")
    print(result_df['grade'].value_counts().sort_index())
    
    print("\nFirst few rows of processed data:")
    print(result_df.head().to_string())
    
    print("\nSample of different stage mappings:")
    stage_mapping = pd.DataFrame({
        'original_stage': df['tumor_stage'].dropna().unique()
    })
    stage_mapping['converted_grade'] = stage_mapping['original_stage'].apply(convert_stage_to_grade)
    print(stage_mapping.sort_values('converted_grade').to_string())
    
if __name__ == "__main__":
    input_file = "/Users/stanleychen/git/Melanoma/data/filtered_patient_data_all.csv"
    output_file = "filtered_patient_data_all_v2.csv"
    
    try:
        process_patient_data(input_file, output_file)
    except Exception as e:
        print(f"An error occurred: {str(e)}")