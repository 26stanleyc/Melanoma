import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

def analyze_patient_data(input_path, output_path):
    """
    Read and analyze patient data with OBS calculation and survival time filtering.
    """
    try:
        # Read all lines first
        with open(input_path, 'r') as file:
            lines = file.readlines()
            
        # Get the first line up to the second occurrence of bcr_patient_uuid
        first_header = lines[0].split('bcr_patient_uuid')[0] + 'bcr_patient_uuid' + lines[0].split('bcr_patient_uuid')[1]
        headers = first_header.strip().split('\t')
        
        # Make headers unique
        unique_headers = []
        seen_counts = {}
        for header in headers:
            if header in seen_counts:
                seen_counts[header] += 1
                unique_headers.append(f"{header}_{seen_counts[header]}")
            else:
                seen_counts[header] = 0
                unique_headers.append(header)
        
        # Read the data with unique headers
        data = pd.read_csv(input_path,
                          sep='\t',
                          names=unique_headers,
                          skiprows=[1, 2])
        
        results = []
        for _, row in data.iterrows():
            try:
                barcode = row['bcr_patient_barcode']
                
                # Get last contact days (OBS)
                try:
                    obs_days = pd.to_numeric(row['last_contact_days_to'], errors='coerce')
                except (ValueError, TypeError):
                    obs_days = None
                
                # Get death days for event status
                try:
                    death_days = pd.to_numeric(row['death_days_to'], errors='coerce')
                    # If death_days exists and is less than or equal to last contact days,
                    # use death_days as the observation time
                    if pd.notna(death_days):
                        if pd.isna(obs_days) or death_days <= obs_days:
                            obs_days = death_days
                except (ValueError, TypeError):
                    death_days = None
                
                # Skip if no valid survival time or negative survival time
                if pd.isna(obs_days) or obs_days < 0:
                    continue
                
                # Find the stage information
                stage_columns = ['stage', 'pathologic_stage', 'clinical_stage', 'Stage', 'ajcc_pathologic_tumor_stage']
                tumor_stage = '[Not Available]'
                
                for col in stage_columns:
                    if col in row.index and pd.notna(row[col]):
                        stage_value = str(row[col]).strip()
                        if stage_value not in ['[Not Available]', '[Not Applicable]', '']:
                            tumor_stage = stage_value
                            break
                
                # Store results
                results.append({
                    'patient_barcode': barcode,
                    'obs_days': obs_days,
                    'tumor_stage': tumor_stage,
                    'vital_status': row['vital_status'],
                    'event': 1 if row['vital_status'] == 'Dead' else 0
                })
                
            except Exception as e:
                print(f"Error processing row for patient {barcode}: {e}")
                continue
        
        # Convert to DataFrame and sort by barcode
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df = results_df.sort_values('patient_barcode')
            
            # Additional filtering to ensure no negative or missing values
            initial_count = len(results_df)
            results_df = results_df.dropna(subset=['obs_days'])
            results_df = results_df[results_df['obs_days'] > 0]
            
            print(f"\nFiltering Results:")
            print(f"Initial patients: {initial_count}")
            print(f"Patients after removing missing/negative survival times: {len(results_df)}")
            print(f"Removed {initial_count - len(results_df)} patients")
            
            # Save filtered results
            results_df.to_csv(output_path, index=False)
            print(f"\nSaved filtered results to: {output_path}")
            
            # Print summary statistics
            print("\nAnalysis Results:")
            print(f"Total patients in final dataset: {len(results_df)}")
            
            print("\nTumor Stage Distribution:")
            print(results_df['tumor_stage'].value_counts())
            
            print("\nVital Status Distribution:")
            print(results_df['vital_status'].value_counts())
            
            valid_obs = results_df['obs_days']
            print(f"\nObserved Survival Interval (OBS) Statistics (in days):")
            print(f"Mean: {valid_obs.mean():.2f}")
            print(f"Median: {valid_obs.median():.2f}")
            print(f"Min: {valid_obs.min():.2f}")
            print(f"Max: {valid_obs.max():.2f}")
            
            event_rate = results_df['event'].mean() * 100
            print(f"\nEvent rate: {event_rate:.1f}%")
        
        return results_df
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        return None

# Define paths
input_path = '/Users/stanleychen/git/Melanoma/clinical_data/nationwidechildrens.org_clinical_patient_skcm.txt'
output_path = '/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv'

# Run analysis
results = analyze_patient_data(input_path, output_path)