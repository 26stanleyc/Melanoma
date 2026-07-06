import pandas as pd
import os
import glob

# Set the path to your main directory
main_dir = '/Users/stanleychen/git/Melanoma/data/combined_data'

# Initialize dictionaries to store data and tracking information
data = {}
missing_data_tracker = {}

# First pass: collect all possible peptide targets
all_peptide_targets = set()
for patient_dir in os.listdir(main_dir):
    patient_path = os.path.join(main_dir, patient_dir)
    if os.path.isdir(patient_path):
        rna_seq_file = glob.glob(os.path.join(patient_path, '*_RNA-seq.tsv'))
        if rna_seq_file:
            with open(rna_seq_file[0], 'r') as f:
                lines = f.readlines()[1:]
            for line in lines:
                peptide_target, _ = line.strip().split(',')
                if peptide_target:
                    all_peptide_targets.add(peptide_target)

# Initialize missing data tracker for each peptide target
for target in all_peptide_targets:
    missing_data_tracker[target] = []

# Second pass: collect data and track missing values
for patient_dir in os.listdir(main_dir):
    patient_path = os.path.join(main_dir, patient_dir)
    if os.path.isdir(patient_path):
        rna_seq_file = glob.glob(os.path.join(patient_path, '*_RNA-seq.tsv'))
        if rna_seq_file:
            # Extract patient ID from the directory name instead of filename
            patient_id = patient_dir  # Use the directory name as the patient ID
            
            # Keep track of which targets we found in this file
            found_targets = set()
            
            # Read the TSV file
            with open(rna_seq_file[0], 'r') as f:
                lines = f.readlines()[1:]
            
            # Process each line
            for line in lines:
                peptide_target, value = line.strip().split(',')
                if peptide_target:
                    found_targets.add(peptide_target)
                    if peptide_target not in data:
                        data[peptide_target] = {}
                    data[peptide_target][patient_id] = float(value)
            
            # Check which targets were missing in this file
            missing_targets = all_peptide_targets - found_targets
            for target in missing_targets:
                missing_data_tracker[target].append(patient_dir)

# Create a DataFrame from the dictionary
result_df = pd.DataFrame.from_dict(data, orient='index')

# Add diagnostic prints
print("\nDiagnostic Information:")
print("-" * 50)
print(f"Number of directories found: {len([d for d in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, d))])}")
all_files = []
for patient_dir in os.listdir(main_dir):
    patient_path = os.path.join(main_dir, patient_dir)
    if os.path.isdir(patient_path):
        rna_seq_file = glob.glob(os.path.join(patient_path, '*_RNA-seq.tsv'))
        if rna_seq_file:
            all_files.extend(rna_seq_file)
print(f"Number of RNA-seq files found: {len(all_files)}")
print(f"Number of unique patient IDs: {len(set(os.path.basename(os.path.dirname(f)) for f in all_files))}")

# Print overall statistics
print("\nOverall Statistics:")
print("-" * 50)
print(f"Total number of peptide targets: {len(all_peptide_targets)}")
print(f"Total number of samples: {len(result_df.columns)}")
print(f"Number of targets with missing values: {sum(1 for x in missing_data_tracker.values() if x)}")

# Save the DataFrame to a CSV file
result_df.to_csv('combined_RNA_data.csv')