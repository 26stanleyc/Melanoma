def merge_expression_files(file1_path, file2_path, output_path):
    # Read first file (log2 TPM values)
    gene_values = {}
    with open(file1_path) as f:
        next(f)  # Skip header
        for _ in range(4):  # Skip first 4 empty rows
            next(f)
        for line in f:
            if line.strip():
                gene, value = line.strip().split(',')
                gene_values[gene] = float(value)
    
    # Read second file and merge
    with open(file2_path) as f, open(output_path, 'w') as out:
        header = f.readline().strip()
        out.write(f"{header}\tTCGA_SKCM\n")
        
        for line in f:
            fields = line.strip().split('\t')
            gene = fields[0]
            if gene in gene_values:  # Only keep overlapping genes
                out.write(f"{line.strip()}\t{gene_values[gene]}\n")

# Example usage:
TCGA_file='/Users/stanleychen/git/Melanoma/final_data/control_folder/control_sample_processed.tsv'
GTEx_file='/Users/stanleychen/git/Melanoma/final_data/control_folder/control_sample_after_gene_id_to_name.tsv'
merge_expression_files(TCGA_file, GTEx_file, 'merged_expression.tsv')