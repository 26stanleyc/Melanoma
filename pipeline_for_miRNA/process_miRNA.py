import os
import pandas as pd

# Directories for input and output
root_directory = "/Users/stanleychen/git/Melanoma/data/miRNA"
processed_directory = "/Users/stanleychen/git/Melanoma/data/miRNA_processed"

# Create the processed directory if it doesn't exist
os.makedirs(processed_directory, exist_ok=True)

# Function to process files recursively
def process_files_recursively(current_directory):
    for entry in os.listdir(current_directory):
        entry_path = os.path.join(current_directory, entry)
        if os.path.isdir(entry_path):  # If entry is a directory, recurse
            process_files_recursively(entry_path)
        elif entry.endswith(".txt"):  # Process only .txt files
            try:
                # Attempt to read the file with different delimiters
                for delimiter in ["\t", " ", ","]:  # Try tab, space, and comma
                    try:
                        df = pd.read_csv(entry_path, delimiter=delimiter)
                        break  # If successful, break out of the loop
                    except pd.errors.ParserError:
                        continue  # Try the next delimiter
                else:
                    print(f"Skipping {entry_path}: Unable to parse with supported delimiters.")
                    continue
                
                # Extract only `miRNA_ID` and `reads_per_million_miRNA_mapped`
                if 'miRNA_ID' in df.columns and 'reads_per_million_miRNA_mapped' in df.columns:
                    processed_df = df[['miRNA_ID', 'reads_per_million_miRNA_mapped']]
                    
                    # Save the processed file directly in the processed directory
                    output_filename = os.path.basename(entry_path)  # Extract file name only
                    processed_path = os.path.join(processed_directory, output_filename)
                    
                    # Save the processed file
                    processed_df.to_csv(processed_path, index=False)
                    print(f"Processed file saved: {processed_path}")
                else:
                    print(f"Skipping {entry_path} as it does not have the required columns.")
            except Exception as e:
                print(f"Error processing {entry_path}: {e}")

# Start processing from the root directory
process_files_recursively(root_directory)
