def count_overlapping_genes(file1_path, file2_path):
    # Read and extract genes from first file
    genes1 = set()
    with open(file1_path, 'r') as f:
        for line in f:
            if line.strip():
                genes1.add(line.split('\t')[0])
    genes1.discard('gene')  # Remove header if present
            
    # Read and extract genes from second file
    genes2 = set()
    with open(file2_path, 'r') as f:
        for line in f:
            if line.strip() and ',' in line:
                gene = line.split(',')[0]
                if gene:  # Only add non-empty gene names
                    genes2.add(gene)
    genes2.discard('')  # Remove empty gene names
    
    # Find overlap
    overlapping_genes = genes1.intersection(genes2)
    return len(overlapping_genes)

# Example usage:
file1='/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/filtered_expression_gene_name.tsv'
file2='/Users/stanleychen/git/Melanoma/final_data/control_folder_v1/control_sample_processed.tsv'
overlap_count = count_overlapping_genes(file1, file2)
print(f"Number of overlapping genes: {overlap_count}")