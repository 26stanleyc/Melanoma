import pandas as pd
import numpy as np
from pathlib import Path

def analyze_patient_data(input_path):
    """
    Read and analyze patient data with proper header handling.
    
    Parameters:
    input_path (str): Path to the clinical data file
    
    Returns:
    pd.DataFrame: DataFrame with patient barcodes, survival times, and tumor status
    """
    try:
        # Read all lines first
        with open(input_path, 'r') as file:
            lines = file.readlines()
            
        # Get the first line up to the second occurrence of bcr_patient_uuid
        first_header = lines[0].split('bcr_patient_uuid')[0] + 'bcr_patient_uuid' + lines[0].split('bcr_patient_uuid')[1]
        headers = first_header.strip().split('\t')
        
        # Read the data, skipping the extra header rows
        data = pd.read_csv(input_path,
                          sep='\t',
                          names=headers,
                          skiprows=[1, 2])  # Skip the second and third header rows
        
        # Create empty lists to store results
        results = []
        
        for _, row in data.iterrows():
            try:
                # Get patient barcode
                barcode = row['bcr_patient_barcode']
                
                # Calculate survival time
                survival_time = None
                survival_source = None
                
                # Check for death_days_to and last_contact_days_to
                if pd.notna(row['death_days_to']):
                    survival_time = row['death_days_to']
                    survival_source = 'death'
                elif pd.notna(row['last_contact_days_to']):
                    survival_time = row['last_contact_days_to']
                    survival_source = 'last_followup'
                
                # Get tumor status
                tumor_status = row['tumor_status']
                if pd.isna(tumor_status):
                    tumor_status = '[Not Available]'
                
                # Store results
                results.append({
                    'patient_barcode': barcode,
                    'survival_time_days': survival_time,
                    'survival_time_source': survival_source,
                    'tumor_status': tumor_status,
                    'vital_status': row['vital_status']
                })
                
            except Exception as e:
                print(f"Error processing row: {e}")
                continue
        
        # Convert to DataFrame and sort by barcode
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df = results_df.sort_values('patient_barcode')
            
            # Print the first few rows for verification
            print("\nFirst few rows of processed data:")
            print(results_df.head())
        
        return results_df
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        # Print available columns for debugging
        if 'data' in locals():
            print("\nAvailable columns:")
            print(data.columns.tolist())
        return None

# Define input path and run analysis
input_data = 'clinical_data.txt'  # Replace with your actual path
results = analyze_patient_data(input_data)

if results is not None and not results.empty:
    # Save results to CSV
    results.to_csv('patient_survival_analysis.csv', index=False)
    
    print("\nAnalysis Results:")
    print(f"Total patients processed: {len(results)}")
    
    print("\nTumor Status Distribution:")
    print(results['tumor_status'].value_counts())
    
    print("\nVital Status Distribution:")
    print(results['vital_status'].value_counts())
    
    print("\nSurvival Time Source Distribution:")
    print(results['survival_time_source'].value_counts())
    
    # Calculate mean survival time only for non-null values
    valid_survival_times = results['survival_time_days'].dropna()
    if len(valid_survival_times) > 0:
        print(f"\nSurvival Time Statistics (in days):")
        print(f"Mean: {valid_survival_times.mean():.2f}")
        print(f"Median: {valid_survival_times.median():.2f}")
        print(f"Min: {valid_survival_times.min():.2f}")
        print(f"Max: {valid_survival_times.max():.2f}")
    else:
        print("\nNo valid survival time data available")
        
    # Save detailed results
    print("\nResults have been saved to 'patient_survival_analysis.csv'")