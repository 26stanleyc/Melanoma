import csv

def get_genes_from_reference(reference_file):
    """Extract gene names from reference file."""
    genes = set()
    with open(reference_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            if row:  # Check if row is not empty
                genes.add(row[0].strip())
    return genes

def filter_patient_data(input_file, output_file, reference_genes):
    """Filter patient data to only include specified genes."""
    with open(input_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        
        # Read and write header
        header = next(reader)
        writer.writerow(header)
        
        # Get gene names from header (first row)
        gene_names = header[1:]  # Skip patient_barcode column
        
        # Create a mask of which columns to keep
        keep_cols = [True]  # Always keep first column (patient_barcode)
        for gene in gene_names:
            # Remove ENST... part if present
            base_gene = gene.split('_')[0]
            keep_cols.append(base_gene in reference_genes)
        
        # Filter and write data
        for row in reader:
            if row:  # Check if row is not empty
                filtered_row = [val for i, val in enumerate(row) if keep_cols[i]]
                writer.writerow(filtered_row)

def main():
    reference_file = "/Users/stanleychen/git/Melanoma/maf_pipeline/top_500_unique_cosmic.csv"  # File with gene list (second image)
    input_file = "/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix_with_0s.csv"        # File with patient data (first image)
    output_file = "/Users/stanleychen/git/Melanoma/maf_pipeline/patient_maf_matrix_500.csv"
    
    # Get reference genes
    reference_genes = get_genes_from_reference(reference_file)
    
    # Filter the data
    filter_patient_data(input_file, output_file, reference_genes)
    print(f"Filtering complete. Results written to {output_file}")

if __name__ == "__main__":
    main()