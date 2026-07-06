import pandas as pd
import numpy as np

def process_patient_data(file_path):
    """
    Process patient data with corrected survival time handling
    """
    # Read the text file with tab delimiter
    df = pd.read_csv(file_path, sep='\t')
    
    # Filter out header-like rows more strictly
    mask = (df['bcr_patient_barcode'].str.contains('^TCGA-', na=False))
    df = df[mask]
    
    # Calculate survival time - Modified logic
    df['survival_time_days'] = None
    
    # For deceased patients, use death_days_to
    deceased_mask = (df['vital_status'] == 'Dead')
    df.loc[deceased_mask, 'survival_time_days'] = df.loc[deceased_mask, 'death_days_to']
    
    # For alive patients, use last_contact_days_to
    alive_mask = (df['vital_status'] == 'Alive')
    df.loc[alive_mask, 'survival_time_days'] = df.loc[alive_mask, 'last_contact_days_to']
    
    # Create source column
    df['survival_time_source'] = np.where(df['vital_status'] == 'Dead', 'death', 'last_followup')
    
    # Debug info
    print("\nBefore final selection:")
    print("Survival times for Dead patients:", df[deceased_mask]['survival_time_days'].notna().sum())
    print("Survival times for Alive patients:", df[alive_mask]['survival_time_days'].notna().sum())
    
    # Select columns for final output
    final_df = pd.DataFrame({
        'patient_barcode': df['bcr_patient_barcode'],
        'survival_time_days': df['survival_time_days'],
        'survival_time_source': df['survival_time_source'],
        'tumor_stage': df['ajcc_pathologic_tumor_stage'],
        'vital_status': df['vital_status']
    })
    
    # Clean up the data
    final_df = final_df.dropna(subset=['survival_time_days'])
    
    # Convert survival time to numeric and round to 1 decimal place
    final_df['survival_time_days'] = pd.to_numeric(final_df['survival_time_days'], errors='coerce')
    final_df = final_df.sort_values('survival_time_days', ascending=False)
    final_df['survival_time_days'] = final_df['survival_time_days'].round(1)
    
    # Additional check for missing survival times
    print("\nFinal checks:")
    print("Missing survival times for Dead patients:", 
          len(final_df[(final_df['vital_status'] == 'Dead') & (final_df['survival_time_days'].isna())]))
    print("Missing survival times for Alive patients:", 
          len(final_df[(final_df['vital_status'] == 'Alive') & (final_df['survival_time_days'].isna())]))
    
    return final_df

def save_filtered_data(df, output_path):
    """
    Save the filtered data to a CSV file
    """
    df.to_csv(output_path, index=False)
    print(f"\nSaved {len(df)} records to {output_path}")

def validate_patient(df, patient_id):
    """
    Validate data for a specific patient
    """
    if patient_id in df['patient_barcode'].values:
        patient_data = df[df['patient_barcode'] == patient_id].iloc[0]
        print(f"\nValidation for patient {patient_id}:")
        print(f"Survival time: {patient_data['survival_time_days']}")
        print(f"Source: {patient_data['survival_time_source']}")
        print(f"Status: {patient_data['vital_status']}")
        return True
    return False

if __name__ == "__main__":
    input_file = "/Users/stanleychen/git/Melanoma/clinical_data/nationwidechildrens.org_clinical_patient_skcm.txt"
    output_file = "filtered_patient_data_all.csv"
    
    try:
        # Process the data
        filtered_df = process_patient_data(input_file)
        
        # Save to CSV
        save_filtered_data(filtered_df, output_file)
        
        # Display summary statistics
        print("\nSummary statistics:")
        print(f"Total patients: {len(filtered_df)}")
        print(f"Alive patients: {len(filtered_df[filtered_df['vital_status'] == 'Alive'])}")
        print(f"Deceased patients: {len(filtered_df[filtered_df['vital_status'] == 'Dead'])}")
        
        # Validate specific patient
        validate_patient(filtered_df, 'TCGA-3N-A9WC')
        
    except FileNotFoundError:
        print(f"Error: Could not find input file '{input_file}'")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        raise