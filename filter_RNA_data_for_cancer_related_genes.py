import pandas as pd

def read_rna_seq_data(file_path):
    """Read RNA-seq data from the first file"""
    data = pd.read_csv(file_path, index_col=0)
    return data

def get_melanoma_and_all_genes(census_file):
    """Extract genes associated with melanoma (in any column) or ALL from the census file"""
    df = pd.read_csv(census_file)
    
    # Convert all columns to string type for safe searching
    df = df.astype(str)
    
    # Create a mask for rows containing 'melanoma' in any column or 'ALL' in any column
    melanoma_mask = df.apply(lambda x: x.str.contains('melanoma', case=False, na=False)).any(axis=1)
    all_mask = df.apply(lambda x: x.str.contains(r'\bALL\b', case=True, na=False)).any(axis=1)
    
    # Combine masks with OR operation
    combined_mask = melanoma_mask | all_mask
    
    # Get unique gene symbols that match either criteria
    target_genes = df[combined_mask]['Gene Symbol'].unique().tolist()
    
    return target_genes

def get_mapped_genes(suppl_table, gencode_file):
    """Get genes from supplementary table and map them using GENCODE data"""
    # Read supplementary table
    suppl_df = pd.read_csv(suppl_table, sep='\t')
    
    # Read GENCODE mapping file
    gencode_df = pd.read_csv(gencode_file, sep='\s+', comment='#')
    
    # Convert Ensembl IDs to lowercase for consistent comparison and remove version numbers
    gencode_df['gene_id'] = gencode_df['gene_id'].str.lower().str.split('.').str[0]
    suppl_df['Ensembl ID'] = suppl_df['Ensembl ID'].str.lower().str.split('.').str[0]
    
    # Create mapping dictionary from Ensembl ID to gene name (keeping original case of gene names)
    gene_mapping = dict(zip(gencode_df['gene_id'], gencode_df['gene_name']))
    
    # Map Ensembl IDs to gene names
    mapped_genes = []
    for ensembl_id in suppl_df['Ensembl ID']:
        if ensembl_id in gene_mapping:
            mapped_genes.append(gene_mapping[ensembl_id])
        
    print(f"\nTotal number of genes mapped: {len(mapped_genes)}")
    print("First 10 mapped genes:", mapped_genes[:10])
    
    return mapped_genes
    
    return mapped_genes

def filter_and_save_data(rna_seq_data, gene_list, output_file):
    """Filter RNA-seq data based on gene list and save to file"""
    # Convert gene list to set for faster lookup
    gene_set = set(gene_list)
    
    # Filter data
    filtered_data = rna_seq_data[rna_seq_data.index.isin(gene_set)]
    
    # Save filtered data
    filtered_data.to_csv(output_file)
    return filtered_data

def main():
    # File paths
    rna_seq_file = '/Users/stanleychen/git/Melanoma/final_data/cleaned_RNA-seq_data.csv'
    census_file = '/Users/stanleychen/git/Melanoma/data/cancer_associated_genes/Census_allMon Jan 13 04_25_31 2025.csv'
    suppl_table_file = '/Users/stanleychen/git/Melanoma/data/cancer_associated_genes/suppl_table1.tsv'
    gencode_file = '/Users/stanleychen/git/Melanoma/data/RNA-seq/0aaec814-e30d-489b-928d-62776b8028cd/dd28c353-13f2-4e6a-a2d7-df94352e043f.rna_seq.augmented_star_gene_counts.tsv'
    output_file = 'RNA_cancer_related_genes.csv'
    
    # Read RNA-seq data
    rna_seq_data = read_rna_seq_data(rna_seq_file)
    
    # Get melanoma and ALL genes
    target_genes = get_melanoma_and_all_genes(census_file)
    print(f"Genes found with melanoma or ALL: {len(target_genes)}")
    print("Sample of target genes:", target_genes[:10])
    
    # Get mapped genes from supplementary table
    mapped_genes = get_mapped_genes(suppl_table_file, gencode_file)
    print(f"Mapped genes from supplementary table: {len(mapped_genes)}")
    print("Sample of mapped genes:", mapped_genes[:10])
    
    # Combine gene lists
    all_genes_to_keep = list(set(target_genes + mapped_genes))
    
    # Filter and save data
    filtered_data = filter_and_save_data(rna_seq_data, all_genes_to_keep, output_file)
    
    # Print summary
    print(f"\nSummary:")
    print(f"Original number of genes: {len(rna_seq_data)}")
    print(f"Number of target genes (melanoma or ALL): {len(target_genes)}")
    print(f"Number of mapped genes from supplementary table: {len(mapped_genes)}")
    print(f"Final number of filtered genes: {len(filtered_data)}")

if __name__ == "__main__":
    main()