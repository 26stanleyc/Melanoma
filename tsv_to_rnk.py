import argparse
import csv
import os

def convert_tsv_to_rnk(input_file, output_file=None):
    """
    Convert a TSV file to RNK format.
    
    Args:
        input_file (str): Path to input TSV file
        output_file (str): Optional path to output RNK file. If not provided,
                          will use input filename with .rnk extension
    """
    if output_file is None:
        # Create output filename by replacing .tsv with .rnk
        output_file = os.path.splitext(input_file)[0] + '.rnk'
    
    try:
        with open(input_file, 'r', encoding='utf-8') as tsv_file:
            # Read TSV file
            tsv_reader = csv.reader(tsv_file, delimiter='\t')
            
            # Get header and data
            header = next(tsv_reader)
            data = list(tsv_reader)
            
            # Write RNK file
            with open(output_file, 'w', encoding='utf-8') as rnk_file:
                # Write header
                rnk_file.write('#1.2\n')  # RNK version header
                rnk_file.write(f'{len(data)}\n')  # Number of entries
                
                # Write column names
                rnk_file.write(f'{len(header)}\n')  # Number of columns
                for column in header:
                    rnk_file.write(f'{column}\n')
                
                # Write data
                for row in data:
                    rnk_file.write('\t'.join(row) + '\n')
                
        print(f"Successfully converted {input_file} to {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found")
    except Exception as e:
        print(f"Error during conversion: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Convert TSV file to RNK format')
    parser.add_argument('input_file', help='Input TSV file path')
    parser.add_argument('-o', '--output', help='Output RNK file path (optional)')
    
    args = parser.parse_args()
    convert_tsv_to_rnk(args.input_file, args.output)

if __name__ == '__main__':
    main()