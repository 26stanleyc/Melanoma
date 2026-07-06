def filter_files_to_overlap(file1_path, file2_path, output1_path, output2_path):
    # First pass - collect genes and their line data from file 2
    genes_file2 = {}  # Will store first occurrence of each gene
    
    # Get first occurrence of each gene from file 2
    with open(file2_path) as f:
        header2 = f.readline()
        for line in f:
            gene = line.strip().split('\t')[0].strip()
            if gene and gene not in genes_file2:  # Only keep first occurrence
                genes_file2[gene] = line

    # Collect genes from file 1 in order
    genes_file1_ordered = []  # Keep track of gene order
    genes_file1_set = set()   # For quick lookup
    
    with open(file1_path) as f:
        header1 = f.readline()
        for line in f:
            gene = line.strip().split(',')[0].strip()
            if gene and gene not in genes_file1_set:
                genes_file1_ordered.append(gene)
                genes_file1_set.add(gene)

    # Find intersection while maintaining file1's order
    overlap_genes_ordered = [gene for gene in genes_file1_ordered if gene in genes_file2]
    print(f"Number of overlapping genes: {len(overlap_genes_ordered)}")
    
    # Write filtered file 1
    count1 = 0
    with open(file1_path) as f, open(output1_path, 'w') as out:
        out.write(header1)
        next(f)  # Skip header again
        for line in f:
            gene = line.strip().split(',')[0].strip()
            if gene in genes_file2:  # if in intersection
                out.write(line)
                count1 += 1
    
    # Write filtered file 2 (matching order of file 1)
    count2 = 0
    with open(output2_path, 'w') as out:
        out.write(header2)
        for gene in overlap_genes_ordered:
            out.write(genes_file2[gene])
            count2 += 1
    
    print(f"\nLines written to output file 1: {count1}")
    print(f"Lines written to output file 2: {count2}")
    
    if count1 != count2:
        print("\nWARNING: Output files have different numbers of lines!")
    
    return overlap_genes_ordered

# Example usage:
tumor_data='/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/combined_rna_seq_data.csv'
control_data='/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/filtered_expression_gene_name.tsv'
filter_files_to_overlap(
    tumor_data,
    control_data,
    'filtered_rna_seq-R.csv',
    'filtered_control-R.tsv'
)
# Example usage:
