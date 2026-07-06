import pandas as pd
import numpy as np

# Load GTEx GCT file properly
control_df = pd.read_csv(
    "/Users/stanleychen/git/Melanoma/GSEA/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct",
    sep="\t",
    skiprows=2,  # Skip metadata lines
    index_col=0  # Set gene IDs as index
)

# Check the number of patients (columns excluding the 'Description' column)
num_patients = control_df.shape[1] - 1  # Exclude gene_name/description column

print(f"Total available patients: {num_patients}")

# Randomly select 400 patient samples (excluding the 'Description' column)
if num_patients > 400:
    selected_patients = np.random.choice(control_df.columns[1:], 400, replace=False)  # Exclude 'Description'
    control_df = control_df[["Description"] + selected_patients.tolist()]  # Keep 'Description' + selected patients
else:
    print("Warning: Less than 400 patients available, using all available samples.")

# Display the filtered control dataset
print(control_df.shape)  # Should be (genes, 401) (one extra column for 'Description')
print(control_df.head())

# Save the filtered dataset
control_df.to_csv("filtered_control_data.csv", sep="\t")
