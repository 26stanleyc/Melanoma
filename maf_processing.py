import os
import shutil

def process_somatic_mutation_directory(source_dir, target_dir):
    # Ensure target directory exists
    os.makedirs(target_dir, exist_ok=True)

    # Walk through the source directory
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith('.maf'):
                # Construct the full file path
                file_path = os.path.join(root, file)
                
                # Construct the target file path
                target_file_path = os.path.join(target_dir, file)
                
                # Copy the file to the target directory
                shutil.copy2(file_path, target_file_path)
                print(f"Copied: {file}")

# Set your source and target directories
source_directory = "/Users/stanleychen/git/Melanoma/data/somatic_mutation"
target_directory = "/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed"

# Run the processing function
process_somatic_mutation_directory(source_directory, target_directory)

print("Processing complete.")