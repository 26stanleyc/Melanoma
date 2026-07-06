import csv
import os

def txt_to_csv(input_txt, output_dir='converted_csv'):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the input filename without path
    input_filename = os.path.basename(input_txt)
    # Create output filename by replacing .txt with .csv
    output_filename = os.path.splitext(input_filename)[0] + '.csv'
    # Create full output path
    output_path = os.path.join(output_dir, output_filename)
    
    print(f"Reading from: {input_txt}")
    print(f"Will save to: {output_path}")
    
    try:
        # Read the text file and convert to CSV
        with open(input_txt, 'r') as txt_file:
            # Read content and split by lines
            lines = txt_file.readlines()
            
        # Write to CSV
        with open(output_path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            for line in lines:
                # Split each line by space or tab and write to CSV
                writer.writerow(line.strip().split())
        
        print(f"\nSuccess! File saved to: {output_path}")
        
    except Exception as e:
        print(f"Error: {str(e)}")

# Example usage
if __name__ == "__main__":
    input_file = "/Users/stanleychen/git/Melanoma/R_scripts/GSEA_with_multiple_controls/h.all.v2024.1.Hs.symbols.gmt"  # Replace with your text file name
    txt_to_csv(input_file)