import gzip
import os
from pathlib import Path
import shutil

def decompress_gz_file(gz_path, output_path=None):
    """
    Decompress a .gz file to the specified output path or same directory.
    
    Args:
        gz_path (str): Path to the .gz file
        output_path (str, optional): Path where to save the decompressed file.
                                   If None, saves to same directory.
    """
    gz_path = Path(gz_path)
    
    # If no output path specified, use the same directory and remove .gz extension
    if output_path is None:
        output_path = gz_path.with_suffix('')
    
    try:
        with gzip.open(gz_path, 'rb') as gz_file:
            with open(output_path, 'wb') as out_file:
                shutil.copyfileobj(gz_file, out_file)
        print(f"Successfully decompressed: {gz_path} -> {output_path}")
    except Exception as e:
        print(f"Error decompressing {gz_path}: {str(e)}")

def decompress_directory(directory_path, recursive=True):
    """
    Decompress all .gz files in a directory.
    
    Args:
        directory_path (str): Path to the directory containing .gz files
        recursive (bool): Whether to process subdirectories
    """
    directory_path = Path(directory_path)
    
    if recursive:
        # Walk through all subdirectories
        for root, _, files in os.walk(directory_path):
            for file in files:
                if file.endswith('.gz'):
                    gz_path = Path(root) / file
                    decompress_gz_file(gz_path)
    else:
        # Only process files in the specified directory
        for file in directory_path.glob('*.gz'):
            decompress_gz_file(file)

if __name__ == "__main__":
    # Example usage:
    
    # To decompress a single file:
    # decompress_gz_file("path/to/your/file.gz")
    
    # To decompress all .gz files in a directory (including subdirectories):
    # decompress_directory("path/to/your/directory")
    
    # To decompress all .gz files in a directory (excluding subdirectories):
    # decompress_directory("path/to/your/directory", recursive=False)
    
    # Example with error handling:
    try:
        # Replace with your file or directory path
        path = "/Users/stanleychen/git/Melanoma/data/maf"
        
        if os.path.isfile(path):
            decompress_gz_file(path)
        elif os.path.isdir(path):
            decompress_directory(path)
        else:
            print(f"Path not found: {path}")
            
    except Exception as e:
        print(f"An error occurred: {str(e)}")