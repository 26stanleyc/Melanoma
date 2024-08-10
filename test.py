import os
import pandas as pd
import glob
from pathlib import Path

def process_file_with_positions(file_path, debug=False):
    # Read the entire RNA-seq file, skipping the first 5 rows and setting column names
    df = pd.read_csv(file_path, sep='\t', skiprows=5, header=None)
    
    if debug:
        print(f"File: {file_path}")
        print(f"Original shape: {df.shape}")
    
    # Extract gene names and strip whitespace
    gene_column = df.iloc[:, 0].str.split(',', expand=True)[0].str.strip().str.upper()
    # print(gene_column)
    # Get positions of genes (indices in the DataFrame)
    positions = df.index.tolist()
    
    return gene_column, positions

# Paths to the files
file1_path = Path('/Users/stanleychen/git/Melanoma/test_combined_data/TCGA-3N-A9WC/TCGA-3N-A9WC_RNA-seq.tsv')
file2_path = Path('/Users/stanleychen/git/Melanoma/test_combined_data/TCGA-3N-A9WD/TCGA-3N-A9WD_RNA-seq.tsv')

# Read the gene columns and positions from both files
genes1_df, positions1 = process_file_with_positions(file1_path, debug=True)
genes2_df, positions2 = process_file_with_positions(file2_path, debug=True)

# print(genes1_df)

# Convert to sets
genes1_set = set(genes1_df.dropna().str.upper())
genes2_set = set(genes2_df.dropna().str.upper())

# Find unmatched genes
unmatched_genes_in_genes1 = genes1_set - genes2_set

# print(unmatched_genes_in_genes1)
unmatched_genes_in_genes2 = genes2_set - genes1_set

# Function to find the positions of genes
def find_gene_positions(gene_series, unmatched_genes, positions):
    gene_to_position = dict(zip(gene_series.str.upper(), positions))
    unmatched_positions = {}
    for gene in unmatched_genes:
        if gene in gene_to_position:
            unmatched_positions[gene] = gene_to_position[gene]
    return unmatched_positions

# Get positions of unmatched genes
unmatched_positions1 = find_gene_positions(genes1_df, unmatched_genes_in_genes1, positions1)
unmatched_positions2 = find_gene_positions(genes2_df, unmatched_genes_in_genes2, positions2)

# Print the first unmatched gene and its positions in both files
def print_matching_unmatched_gene(unmatched_positions1, unmatched_positions2):
    for gene in unmatched_positions1:
        if gene in unmatched_positions2:
            print(f"Unmatched gene '{gene}' found in both files.")
            print(f"Position in file 1: {unmatched_positions1[gene]}")
            print(f"Position in file 2: {unmatched_positions2[gene]}")
            break

print_matching_unmatched_gene(unmatched_positions1, unmatched_positions2)

print(f"Number of unmatched genes: {len(unmatched_genes_in_genes1.union(unmatched_genes_in_genes2))} out of {len(genes1_set.union(genes2_set))}")
