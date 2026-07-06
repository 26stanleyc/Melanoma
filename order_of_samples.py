import os
import pandas as pd

import os

rna_seq_path = 'RNA-seq'
dir_list = sorted([d for d in os.listdir(rna_seq_path) if os.path.isdir(os.path.join(rna_seq_path, d))])
print(dir_list)
# Now dir_list will match the file explorer's order since it's using standard 
# lexicographical sorting which matches the hex string ordering we see

mapping_df = pd.read_csv('/Users/stanleychen/git/Melanoma/clinical_data/data/RNA-seq/gdc_sample_sheet.2024-12-27.tsv', sep='\t')

with open('sample_order.txt', 'w') as f:
    for d in dir_list:
        case_id = mapping_df[mapping_df['File ID'] == d]['Case ID'].iloc[0]
        f.write(f"{case_id}\n")