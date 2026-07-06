def count_dimensions(file_path):
    """
    Count total number of rows and columns in a large CSV file efficiently.
    
    Args:
        file_path (str): Path to the CSV file
    Returns:
        tuple: (num_rows, num_columns)
    """
    num_rows = 0
    num_columns = 0
    
    with open(file_path, 'r') as file:
        # Get number of columns from first line
        first_line = file.readline().strip()
        num_columns = len(first_line.split(','))
        
        # Count first line
        num_rows = 1
        
        # Count remaining lines
        for _ in file:
            num_rows += 1
    
    print(f"Number of rows: {num_rows}")
    print(f"Number of columns: {num_columns}")
    
    return num_rows, num_columns


# Usage
if __name__ == "__main__":
    file_path = "/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA_after_IGCBA.csv"  # Replace with your file path
    try:
        num_cols = count_dimensions(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
    