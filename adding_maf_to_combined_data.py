import os
import pandas as pd
import shutil

def process_maf_files(combined_data_dir, somatic_mutation_dir, sample_sheet_path):
    # Read the sample sheet
    sample_sheet = pd.read_csv(sample_sheet_path, sep='\t')
    # Create a dictionary mapping File Names to Case IDs (using only the first part)
    file_to_case = {row['File Name'].replace('.gz', ''): row['Case ID'].split(',')[0].strip() 
                    for _, row in sample_sheet.iterrows()}
    # print(file_to_case)
    # Get the list of folders in combined_data
    case_folders = [f for f in os.listdir(combined_data_dir) if os.path.isdir(os.path.join(combined_data_dir, f))]
    
    print("Folders in combined_data:")
    for folder in case_folders:
        print(folder)
    print(f"Total number of folders: {len(case_folders)}")
    print("\n" + "="*50 + "\n")

    # Iterate through files in the somatic_mutation directory
    for filename in os.listdir(somatic_mutation_dir):
        if filename.endswith('.maf'):
            file_path = os.path.join(somatic_mutation_dir, filename)
            # print(file_path)
            # Find the corresponding Case ID
            case_id = file_to_case.get(filename)
            # print(case_id)
            if case_id:
                # Find matching folder in combined_data
                matching_folder = next((folder for folder in case_folders if case_id in folder), None)
                
                if matching_folder:
                    # Construct the path to the case folder
                    case_folder = os.path.join(combined_data_dir, matching_folder)
                    
                    # Rename the file using the folder name
                    new_filename = f"{matching_folder}_maf"
                    
                    # Copy and rename the file to the case folder
                    destination = os.path.join(case_folder, new_filename)
                    shutil.copy2(file_path, destination)
                    print(f"Copied and renamed {filename} to {destination}")
                else:
                    print(f"No matching folder found for Case ID: {case_id}")
            else:
                print(f"No matching Case ID found for {filename}")


# Set your directory paths
combined_data_dir = "/Users/stanleychen/git/Melanoma/data/combined_data"
somatic_mutation_dir = "/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed"
sample_sheet_path = "data/somatic_mutation/gdc_sample_sheet.2024-09-28.tsv"

# Run the function
process_maf_files(combined_data_dir, somatic_mutation_dir, sample_sheet_path)

print("Processing complete.")