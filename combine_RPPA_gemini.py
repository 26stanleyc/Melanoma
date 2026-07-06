import pandas as pd
import os
import glob

main_dir = '/Users/stanleychen/git/Melanoma/data/combined_data_no_RPPA'
data = {}

for patient_dir in os.listdir(main_dir):
    patient_path = os.path.join(main_dir, patient_dir)
    if os.path.isdir(patient_path):
        rna_seq_file = glob.glob(os.path.join(patient_path, '*_RNA-seq.tsv'))
        if rna_seq_file:
            patient_id = os.path.basename(rna_seq_file[0]).split('_')[0]
            
            with open(rna_seq_file[0], 'r') as f:
                lines = f.readlines()[1:]
            
            for line in lines:
                parts = line.strip().split(',')
                gene = parts[0]
                value = parts[2]  # Take the third column (log2_tpm_plus_1)
                if gene:
                    if gene not in data:
                        data[gene] = {}
                    data[gene][patient_id] = float(value)

result_df = pd.DataFrame.from_dict(data, orient='index')
print(result_df.head())
result_df.to_csv('combined_RNA-seq_data.csv')