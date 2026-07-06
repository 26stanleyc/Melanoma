import pandas as pd
import numpy as np
from scipy import stats

def calculate_ranked_stats(file1_path, file2_path, output_path):
    # Read CSV files
    df1 = pd.read_csv(file1_path, index_col=0)
    df2 = pd.read_csv(file2_path, index_col=0)
    
    results_list = []
    
    for gene in df1.index:
        if gene in df2.index:
            gene_data1 = df1.loc[gene].dropna()
            gene_data2 = df2.loc[gene].dropna()
            
            if len(gene_data1) > 0 and len(gene_data2) > 0:
                # Calculate t-statistic only
                t_stat, _ = stats.ttest_ind(gene_data1, gene_data2)
                
                results_list.append({
                    'Symbol': gene,
                    'Ranked Stat': t_stat
                })
    
    # Create, sort and round results
    results = pd.DataFrame(results_list)
    results = results.sort_values('Ranked Stat', ascending=False)
    results = results.round(6)
    
    # Save results
    results.to_csv(output_path, sep='\t', index=False)
    
    print(f"\nTotal genes ranked: {len(results)}")
    
    return results

# Example usage
file1_path = "/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/filtered_rna_seq-R.csv"
file2_path = "/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/filtered_expression_log2_tpm-R.csv"
output_path = '/Users/stanleychen/git/Melanoma/R_scripts/GSEA_with_multiple_controls/ranked_list.tsv'

results = calculate_ranked_stats(file1_path, file2_path, output_path)