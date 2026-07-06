import pandas as pd
import numpy as np
import os
import gseapy as gp
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce

# File paths
combined_data_path = "/Users/stanleychen/git/Melanoma/GSEA/combined_case_control.csv"
output_dir = "/Users/stanleychen/git/Melanoma/GSEA/results"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# ===============================
# 1. Load and Process the Data
# ===============================
print("Loading combined dataset...")
df = pd.read_csv(combined_data_path, index_col=0)

# Separate cancer (TCGA) and control (GTEX) columns
cancer_cols = [col for col in df.columns if col.startswith('TCGA')]
control_cols = [col for col in df.columns if col.startswith('GTEX')]

print(f"Found {len(cancer_cols)} cancer samples and {len(control_cols)} control samples")

# ===============================
# 2. Calculate Differential Expression
# ===============================
print("Calculating differential expression...")

# Calculate mean expression for each gene in cancer and control groups
cancer_means = df[cancer_cols].mean(axis=1)
control_means = df[control_cols].mean(axis=1)

# Calculate log2 fold change
log2fc = np.log2((cancer_means + 0.01) / (control_means + 0.01))  # Adding small value to avoid division by zero

# Perform t-test for each gene
pvalues = []
for gene in df.index:
    cancer_vals = df.loc[gene, cancer_cols]
    control_vals = df.loc[gene, control_cols]
    _, pval = stats.ttest_ind(cancer_vals, control_vals, equal_var=False)
    pvalues.append(pval)

# Adjust p-values for multiple testing
padj = stats.false_discovery_control(pvalues)

# Create a DataFrame with differential expression results
diff_exp = pd.DataFrame({
    'log2FoldChange': log2fc,
    'pvalue': pvalues,
    'padj': padj
}, index=df.index)

# Save differential expression results
diff_exp.to_csv(os.path.join(output_dir, 'differential_expression.csv'))

# Filter for significantly differentially expressed genes
sig_genes = diff_exp[(diff_exp['padj'] < 0.05) & (abs(diff_exp['log2FoldChange']) > 1)]
print(f"Found {len(sig_genes)} significantly differentially expressed genes")

# ===============================
# 3. Prepare Pre-ranked Gene List for GSEA
# ===============================
print("Preparing pre-ranked gene list for GSEA...")

# Rank genes by sign of fold change * -log10(pvalue)
diff_exp['score'] = diff_exp['log2FoldChange'] * -np.log10(diff_exp['pvalue'] + 1e-10)
ranked_genes = diff_exp.sort_values('score', ascending=False)

# Save the ranked gene list for GSEA
rnk_file = os.path.join(output_dir, 'gsea_ranked_genes.rnk')
ranked_genes[['score']].to_csv(rnk_file, sep='\t', header=False)

# ===============================
# 4. Run GSEA Pre-ranked Analysis
# ===============================
print("Running GSEA...")

# Collections to use
gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021', 'Reactome_2022']

for gene_set in gene_sets:
    print(f"Running GSEA with {gene_set}...")
    gsea_result = gp.prerank(
        rnk=rnk_file,
        gene_sets=gene_set,
        outdir=os.path.join(output_dir, f'GSEA_{gene_set}'),
        min_size=10,
        max_size=500,
        permutation_num=1000,
        threads=4,  # Updated from processes to threads
        seed=42,
        no_plot=False
    )
    
    # Save the top enriched pathways to a separate file
    top_pathways = gsea_result.res2d.sort_values('NES', ascending=False)
    top_pathways.to_csv(os.path.join(output_dir, f'top_pathways_{gene_set}.csv'))

# ===============================
# 5. Generate Summary Visualizations
# ===============================
print("Generating summary visualizations...")

# 1. Volcano plot of differential expression
plt.figure(figsize=(10, 8))
plt.scatter(
    diff_exp['log2FoldChange'], 
    -np.log10(diff_exp['pvalue']), 
    alpha=0.5, 
    color='gray', 
    s=20
)

# Highlight significant genes
sig_genes_plot = diff_exp[(diff_exp['padj'] < 0.05) & (abs(diff_exp['log2FoldChange']) > 1)]
plt.scatter(
    sig_genes_plot['log2FoldChange'], 
    -np.log10(sig_genes_plot['pvalue']), 
    alpha=0.8, 
    color='red', 
    s=25
)

# Add some gene labels for top genes
top_genes = sig_genes_plot.sort_values('score', ascending=False).head(15)
for gene in top_genes.index:
    plt.annotate(
        gene, 
        xy=(diff_exp.loc[gene, 'log2FoldChange'], -np.log10(diff_exp.loc[gene, 'pvalue'])),
        xytext=(5, 5), 
        textcoords='offset points',
        fontsize=8
    )

plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
plt.axvline(-1, linestyle='--', color='gray')
plt.axvline(1, linestyle='--', color='gray')
plt.xlabel('log2 Fold Change (Cancer / Control)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot of Differential Expression')
plt.savefig(os.path.join(output_dir, 'volcano_plot.png'), dpi=300, bbox_inches='tight')

# 2. Combined enrichment plot showing top pathways across all gene sets
all_top_pathways = []
for gene_set in gene_sets:
    try:
        top_path = pd.read_csv(os.path.join(output_dir, f'top_pathways_{gene_set}.csv'), index_col=0)
        top_path['Gene_Set_Collection'] = gene_set
        top_path_subset = top_path[top_path['FDR q-val'] < 0.25].head(5)  # Top 5 with FDR < 0.25
        all_top_pathways.append(top_path_subset)
    except:
        print(f"Could not load results for {gene_set}. Skipping.")

if all_top_pathways:
    combined_top = pd.concat(all_top_pathways)
    combined_top = combined_top.sort_values('NES', ascending=False)
    
    plt.figure(figsize=(12, 8))
    barplot = sns.barplot(
        x='NES', 
        y=combined_top.index, 
        hue='Gene_Set_Collection',
        data=combined_top,
        palette='viridis'
    )
    plt.title('Top Enriched Pathways (NES Score)')
    plt.xlabel('Normalized Enrichment Score')
    plt.ylabel('Pathway')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'top_pathways_combined.png'), dpi=300, bbox_inches='tight')

# ===============================
# 6. Create Summary Report
# ===============================
print("Creating summary report...")

# Create a simple HTML report
html_report = f"""
<!DOCTYPE html>
<html>
<head>
    <title>GSEA Results: Cancer vs Control</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1, h2, h3 {{ color: #333366; }}
        .results-section {{ margin: 20px 0; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
    </style>
</head>
<body>
    <h1>Gene Set Enrichment Analysis: Cancer vs Control</h1>
    
    <div class="results-section">
        <h2>Differential Expression Summary</h2>
        <p>Total genes analyzed: {len(diff_exp)}</p>
        <p>Significantly differentially expressed genes (padj < 0.05, |log2FC| > 1): {len(sig_genes)}</p>
        <p>Top upregulated genes in Cancer:</p>
        <table>
            <tr><th>Gene</th><th>log2 Fold Change</th><th>Adjusted p-value</th></tr>
            {reduce(lambda x, y: x + y, ['<tr><td>{0}</td><td>{1:.3f}</td><td>{2:.3e}</td></tr>'.format(
                gene, diff_exp.loc[gene, 'log2FoldChange'], diff_exp.loc[gene, 'padj']
            ) for gene in diff_exp.sort_values('log2FoldChange', ascending=False).head(10).index])}
        </table>
        
        <p>Top downregulated genes in Cancer:</p>
        <table>
            <tr><th>Gene</th><th>log2 Fold Change</th><th>Adjusted p-value</th></tr>
            {reduce(lambda x, y: x + y, ['<tr><td>{0}</td><td>{1:.3f}</td><td>{2:.3e}</td></tr>'.format(
                gene, diff_exp.loc[gene, 'log2FoldChange'], diff_exp.loc[gene, 'padj']
            ) for gene in diff_exp.sort_values('log2FoldChange', ascending=True).head(10).index])}
        </table>
    </div>
    
    <div class="results-section">
        <h2>GSEA Results</h2>
        <h3>Top Enriched Pathways</h3>
"""

# Add results for each gene set collection
for gene_set in gene_sets:
    try:
        top_path = pd.read_csv(os.path.join(output_dir, f'top_pathways_{gene_set}.csv'), index_col=0)
        
        html_report += f"""
        <h4>{gene_set}</h4>
        <table>
            <tr><th>Pathway</th><th>NES</th><th>FDR q-val</th><th>Gene Set Size</th></tr>
            {reduce(lambda x, y: x + y, ['<tr><td>{0}</td><td>{1:.3f}</td><td>{2:.3e}</td><td>{3}</td></tr>'.format(
                pathway, row['NES'], row['FDR q-val'], int(row['Size'])
            ) for pathway, row in top_path.sort_values('NES', ascending=False).head(10).iterrows()])}
        </table>
        """
    except:
        html_report += f"<p>No results available for {gene_set}</p>"

html_report += """
    </div>
    
    <div class="results-section">
        <h2>Visualizations</h2>
        <h3>Volcano Plot</h3>
        <img src="volcano_plot.png" alt="Volcano Plot" style="max-width:100%; height:auto;">
        
        <h3>Top Pathways</h3>
        <img src="top_pathways_combined.png" alt="Top Pathways" style="max-width:100%; height:auto;">
    </div>
    
    <div class="results-section">
        <h2>Methods</h2>
        <p>This analysis used the following steps:</p>
        <ol>
            <li>Loaded combined cancer and control expression data</li>
            <li>Calculated differential expression using t-tests with FDR correction</li>
            <li>Ranked genes by sign(log2FC) * -log10(pvalue)</li>
            <li>Performed GSEA using pre-ranked algorithm</li>
            <li>Generated summary visualizations and reports</li>
        </ol>
    </div>
</body>
</html>
"""

with open(os.path.join(output_dir, 'gsea_summary_report.html'), 'w') as f:
    f.write(html_report)

print(f"GSEA analysis complete! Results saved to {output_dir}")
print(f"Check the summary report at {os.path.join(output_dir, 'gsea_summary_report.html')}")