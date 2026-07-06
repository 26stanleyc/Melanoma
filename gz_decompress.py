import os
import gzip
import shutil

def process_folder(root_folder):
    for dirpath, dirnames, filenames in os.walk(root_folder):
        for filename in filenames:
            if filename.endswith('.gz'):
                file_path = os.path.join(dirpath, filename)
                output_path = os.path.join(dirpath, filename[:-3])  # Remove .gz extension
                
                with gzip.open(file_path, 'rb') as f_in:
                    with open(output_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                
                print(f"Decompressed: {file_path}")
            else:
                print(f"Skipped: {os.path.join(dirpath, filename)}")

# Replace with your actual folder path
root_folder = '/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed_v3/'
process_folder(root_folder)