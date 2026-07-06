import os
import pandas as pd
import glob
import shutil

# Read sample sheets
rna_samples = pd.read_csv('/Users/stanleychen/git/Melanoma/data/RNA-seq/gdc_sample_sheet.2024-12-27.tsv', sep='\t')
mirna_samples = pd.read_csv('/Users/stanleychen/git/Melanoma/data/miRNA/gdc_sample_sheet.2025-02-03.tsv', sep='\t')

# Extract Case IDs
rna_case_ids = set(rna_samples['Case ID'])
mirna_case_ids = set(id.split(',')[0].strip() for id in mirna_samples['Case ID'])

# Find matching Case IDs
matching_case_ids = rna_case_ids.intersection(mirna_case_ids)
print(f"Matching Case IDs:", len(matching_case_ids))

def process_files(rna_samples, mirna_samples, case_ids):
    for case_id in case_ids:
        rna_row = rna_samples[rna_samples['Case ID'] == case_id]
        mirna_row = mirna_samples[mirna_samples['Case ID'] == case_id]
        
        if not rna_row.empty and not mirna_row.empty:
            rna_file_name = rna_row['File Name'].values[0]
            mirna_file_name = mirna_row['File Name'].values[0]
            
            rna_file_name = rna_file_name.replace('.tsv', '_processed.tsv')
            mirna_file_name = mirna_file_name.replace('.tsv', '_processed.tsv')
            
            # Find the files
            rna_files = glob.glob(os.path.join('RNA-seq_processed', '**', rna_file_name), recursive=True)
            mirna_files = glob.glob(os.path.join('miRNA_processed', '**', mirna_file_name), recursive=True)
            
            # Check if both files were found
            if rna_files and mirna_files:
                rna_file_path = rna_files[0]  # Take the first match
                mirna_file_path = mirna_files[0]  # Take the first match
                
                # Create case folder
                case_folder = os.path.join('combined_data_no_RPPA', case_id)
                os.makedirs(case_folder, exist_ok=True)
                
                # Copy files
                rna_dest = os.path.join(case_folder, f"{case_id}_RNA-seq.tsv")
                mirna_dest = os.path.join(case_folder, f"{case_id}_miRNA.tsv")
                
                try:
                    shutil.copy2(rna_file_path, rna_dest)
                    shutil.copy2(mirna_file_path, mirna_dest)
                    print(f"Created folder and copied files for Case ID: {case_id}")
                except Exception as e:
                    print(f"Error processing Case ID {case_id}: {str(e)}")

# Process files
process_files(rna_samples, mirna_samples, matching_case_ids)
print("Finished processing and organizing files.")