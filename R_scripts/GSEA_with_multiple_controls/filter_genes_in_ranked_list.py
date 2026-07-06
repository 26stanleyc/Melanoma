import pandas as pd

def read_gmt_file(gmt_file_path):
    """Read GMT file and return set of all unique genes."""
    all_genes = set()
    
    with open(gmt_file_path, 'r') as file:
        for line in file:
            # Split the line into parts
            parts = line.strip().split('\t')
            # Skip pathway name and URL (first two elements)
            genes = parts[2:]
            # Add genes to set
            all_genes.update(genes)
    
    return all_genes

def filter_ranked_genes(ranked_file_path, gmt_file_path, output_path):
    """Filter ranked gene file to only include genes present in GMT file."""
    # Read the ranked gene file
    ranked_df = pd.read_csv(ranked_file_path, sep='\t')
    
    # Get all unique genes from GMT file
    hallmark_genes = read_gmt_file(gmt_file_path)
    
    # Filter ranked file to only include genes in hallmark set
    filtered_df = ranked_df[ranked_df['Symbol'].isin(hallmark_genes)]
    
    # Sort by Ranked Stat (assuming this is the column name)
    filtered_df = filtered_df.sort_values('Ranked Stat', ascending=False)
    
    # Save filtered results
    filtered_df.to_csv(output_path, sep='\t', index=False)
    
    # Print some statistics
    total_genes = len(ranked_df)
    filtered_genes = len(filtered_df)
    print(f"Original number of genes: {total_genes}")
    print(f"Number of genes after filtering: {filtered_genes}")
    print(f"Removed {total_genes - filtered_genes} genes")
    
    return filtered_df

# Example usage
ranked_file = "/Users/stanleychen/git/Melanoma/R_scripts/GSEA_with_multiple_controls/ranked_list.tsv"  # Your ranked gene file
gmt_file = "/Users/stanleychen/git/Melanoma/R_scripts/GSEA_with_multiple_controls/h.all.v2024.1.Hs.symbols.gmt"       # Your hallmark GMT file
output_file = "filtered_ranked_list.tsv"

filtered_results = filter_ranked_genes(ranked_file, gmt_file, output_file)