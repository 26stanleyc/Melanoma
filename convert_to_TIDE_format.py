import pandas as pd
import os
import glob

# Set the path to your main directory
main_dir = '/Users/stanleychen/git/Melanoma/data/combined_data'

# Initialize an empty dictionary to store the data
data = {}

# Walk through the directory structure
for patient_dir in os.listdir(main_dir):
    patient_path = os.path.join(main_dir, patient_dir)
    if os.path.isdir(patient_path):
        # Look for the RNA-seq file
        rna_seq_file = glob.glob(os.path.join(patient_path, '*_RNA-seq.tsv'))
        if rna_seq_file:
            # Extract patient ID from the filename
            patient_id = os.path.basename(rna_seq_file[0]).split('_')[0]
            
            # Read the TSV file
            with open(rna_seq_file[0], 'r') as f:
                lines = f.readlines()[1:]  # Skip the header
            
            # Process each line
            for line in lines:
                gene, value = line.strip().split(',')
                if gene:  # Only process lines with a gene name
                    if gene not in data:
                        data[gene] = {}
                    data[gene][patient_id] = float(value)

# Create a DataFrame from the dictionary
result_df = pd.DataFrame.from_dict(data, orient='index')

# Display the first few rows of the DataFrame
print(result_df.head())

# Save the DataFrame to a CSV file
result_df.to_csv('combined_rna_seq_data.csv')