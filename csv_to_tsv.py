import csv

def csv_to_tsv(csv_file, tsv_file):
    """
    Convert a CSV file to a TSV file.

    Args:
    csv_file (str): Path to the input CSV file
    tsv_file (str): Path to the output TSV file

    Returns:
    None
    """
    try:
        with open(csv_file, 'r', newline='') as csv_input, \
             open(tsv_file, 'w', newline='') as tsv_output:
            
            csv_reader = csv.reader(csv_input)
            tsv_writer = csv.writer(tsv_output, delimiter='\t')
            
            for row in csv_reader:
                tsv_writer.writerow(row)
        
        print(f"Successfully converted {csv_file} to {tsv_file}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

# Example usage
csv_to_tsv('normalized_rna_seq_data.csv', 'normalized_rna_seq_data.tsv')