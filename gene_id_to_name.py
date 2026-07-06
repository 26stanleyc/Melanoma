def create_gene_mapping(annotation_file):
    mapping = {}
    with open(annotation_file, 'r') as f:
        for i, line in enumerate(f):
            if i < 6:  # Skip first 6 rows
                continue
            fields = line.strip().split('\t')
            gene_id = fields[0].split('.')[0]  # Remove version number
            gene_name = fields[1]
            mapping[gene_id] = gene_name
    return mapping

def convert_genes_in_file(input_file, output_file, gene_mapping):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # Write header
        header = next(f_in)
        f_out.write(header)
        
        # Process data lines
        for line in f_in:
            fields = line.strip().split('\t')
            gene_id = fields[0]
            
            # If gene_id starts with ENSG, try to convert it
            if gene_id.startswith('ENSG'):
                base_gene_id = gene_id.split('.')[0]
                if base_gene_id in gene_mapping:
                    fields[0] = gene_mapping[base_gene_id]
            
            f_out.write('\t'.join(fields) + '\n')

# Example usage:
mapping = create_gene_mapping('/Users/stanleychen/git/Melanoma/final_data/control_folder_v1/control_sample.tsv')
convert_genes_in_file('/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/filtered_expression.txt', 'output_data.tsv', mapping)