import csv
from collections import OrderedDict

def process_gene_names(input_file, output_file, limit=500):
    # Use OrderedDict to maintain insertion order while ensuring uniqueness
    unique_genes = OrderedDict()
    
    with open(input_file, 'r') as f:
        # Skip header row
        next(f)
        
        reader = csv.reader(f)
        for row in reader:
            if not row:  # Skip empty rows
                continue
                
            # Get the gene name (first column)
            full_gene_name = row[0].strip()
            
            # Extract base gene name by splitting on underscore and taking first part
            base_gene_name = full_gene_name.split('_')[0]
            
            # Add to ordered dict if we haven't seen it before
            if base_gene_name not in unique_genes and len(unique_genes) < limit:
                unique_genes[base_gene_name] = ','.join(row[1:])  # Keep the other columns
    
    # Write the results to output file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow(['Gene name', 'Mutated samples', 'Samples tested'])
        
        # Write unique gene entries
        for gene_name, other_cols in unique_genes.items():
            writer.writerow([gene_name] + other_cols.split(','))

    print(f"Processed {len(unique_genes)} unique genes")
    return len(unique_genes)

# Example usage
if __name__ == "__main__":
    input_file = "/Users/stanleychen/git/Melanoma/data/Genes_with_mutationTue Feb 11 22_45_05 2025.csv"  # Replace with your input file name
    output_file = "/Users/stanleychen/git/Melanoma/maf_pipeline/top_500_unique_cosmic.csv"  # Replace with desired output file name
    num_processed = process_gene_names(input_file, output_file)