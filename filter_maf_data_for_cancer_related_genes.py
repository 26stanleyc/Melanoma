import pandas as pd

def read_variant_data(file_path):
    """Read the variant data from the input CSV file"""
    # Read the CSV with pandas, treating first row as header
    df = pd.read_csv(file_path)
    
    # Print column names to debug
    print("\nColumns in variant data:")
    print(df.columns.tolist())
    
    return df

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
    
    # Get unique gene symbols and convert to uppercase
    target_genes = [gene.upper() for gene in df[combined_mask]['Gene Symbol'].unique().tolist()]
    
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
    
    # Create mapping dictionary from Ensembl ID to gene name
    gene_mapping = dict(zip(gencode_df['gene_id'], gencode_df['gene_name']))
    
    # Map Ensembl IDs to gene names and convert to uppercase
    mapped_genes = []
    for ensembl_id in suppl_df['Ensembl ID']:
        if ensembl_id in gene_mapping:
            mapped_genes.append(gene_mapping[ensembl_id].upper())
    
    return mapped_genes

def filter_and_save_data(variant_data, gene_list, output_file):
    """Filter variant data based on gene list and save to file"""
    # Print column names and first few rows for debugging
    print("\nFirst few rows of variant data:")
    print(variant_data.head())
    
    # Convert second column to uppercase for comparison (assuming gene is second column)
    variant_data['gene_upper'] = variant_data.iloc[:, 1].str.upper()
    gene_set = set(gene.upper() for gene in gene_list)
    
    # Print some debug information
    print("\nDebugging information:")
    print("Sample of variant genes:", sorted(variant_data.iloc[:, 1].unique())[:10])
    print("Sample of gene list:", sorted(list(gene_set))[:10])
    
    # Check for any matches
    matches = variant_data['gene_upper'].isin(gene_set)
    print(f"Number of matching genes found: {matches.sum()}")
    
    # Filter data
    filtered_data = variant_data[matches].drop(columns=['gene_upper'])
    
    # Save filtered data
    filtered_data.to_csv(output_file, index=False)
    return filtered_data

def main():
    # File paths
    variant_file = '/Users/stanleychen/git/Melanoma/final_data/mutations_with_alphamissense.csv'
    census_file = '/Users/stanleychen/git/Melanoma/data/cancer_associated_genes/Census_allMon Jan 13 04_25_31 2025.csv'
    suppl_table_file = '/Users/stanleychen/git/Melanoma/data/cancer_associated_genes/suppl_table1.tsv'
    gencode_file = '/Users/stanleychen/git/Melanoma/data/RNA-seq/0aaec814-e30d-489b-928d-62776b8028cd/dd28c353-13f2-4e6a-a2d7-df94352e043f.rna_seq.augmented_star_gene_counts.tsv'
    output_file = 'filtered_mutation_data.csv'
    
    # Read variant data
    variant_data = read_variant_data(variant_file)
    print(f"Original number of variants: {len(variant_data)}")
    
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
    filtered_data = filter_and_save_data(variant_data, all_genes_to_keep, output_file)
    
    # Print summary
    print(f"\nSummary:")
    print(f"Original number of variants: {len(variant_data)}")
    print(f"Number of target genes (melanoma or ALL): {len(target_genes)}")
    print(f"Number of mapped genes from supplementary table: {len(mapped_genes)}")
    print(f"Final number of filtered variants: {len(filtered_data)}")
    
    # Print genes that were actually found in the variant data
    if len(filtered_data) > 0:
        print("\nGenes found in variant data:")
        print(sorted(filtered_data.iloc[:, 1].unique()))

if __name__ == "__main__":
    main()

