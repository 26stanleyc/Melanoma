import os
import pandas as pd
import glob
import shutil

# Read sample sheets
rna_samples = pd.read_csv('data/RNA-seq/gdc_sample_sheet.2024-08-05.tsv', sep='\t')
rppa_samples = pd.read_csv('data/RPPA/gdc_sample_sheet.2024-08-05 (1).tsv', sep='\t')

# Extract Case IDs
rna_case_ids = set(rna_samples['Case ID'])
rppa_case_ids = set(rppa_samples['Case ID'])

# Find matching Case IDs
matching_case_ids = rna_case_ids.intersection(rppa_case_ids)

# print(f"Matching Case IDs: {matching_case_ids}")

# Function to find and copy matched files
def process_files(rna_samples, rppa_samples, case_ids):
    for case_id in case_ids:
        rna_row = rna_samples[rna_samples['Case ID'] == case_id]
        rppa_row = rppa_samples[rppa_samples['Case ID'] == case_id]
        if not rna_row.empty and not rppa_row.empty:
            rna_file_name = rna_row['File Name'].values[0]
            rppa_file_name = rppa_row['File Name'].values[0]
            rna_file_name = rna_file_name.replace('.tsv', '_processed.tsv')
            rppa_file_name = rppa_file_name.replace('.tsv', '_processed.tsv')
            rna_file_path = glob.glob(os.path.join('RNA-seq_processed', '**', rna_file_name), recursive=True)
            rppa_file_path = glob.glob(os.path.join('RPPA_processed', '**', rppa_file_name), recursive=True)
            if rna_file_path and rppa_file_path:
                rna_file_path = rna_file_path[0]
                rppa_file_path = rppa_file_path[0]
                # Create a new folder for this Case ID
                case_folder = os.path.join('combined_data', case_id)
                os.makedirs(case_folder, exist_ok=True)

                # Copy RNA-seq file
                rna_dest = os.path.join(case_folder, f"{case_id}_RNA-seq.tsv")
                shutil.copy2(rna_file_path, rna_dest)

                # Copy RPPA file
                rppa_dest = os.path.join(case_folder, f"{case_id}_RPPA.tsv")
                shutil.copy2(rppa_file_path, rppa_dest)

                print(f"Created folder and copied files for Case ID: {case_id}")

# Process files
process_files(rna_samples, rppa_samples, matching_case_ids)

print("Finished processing and organizing files.")
