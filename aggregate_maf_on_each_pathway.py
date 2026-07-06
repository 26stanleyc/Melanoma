import pandas as pd
import numpy as np
from collections import defaultdict

def load_pathway_genes(pathway_file):
    """Load pathway definitions and their genes."""
    pathways = {}
    
    with open(pathway_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:  # Ensure we have pathway name and genes
                pathway_name = parts[0]
                # Get genes (skip first two columns which are pathway name and URL)
                genes = [gene for gene in parts[2:] if gene.strip()]
                pathways[pathway_name] = set(genes)
    
    return pathways

def load_pathogenicity_scores(scores_file):
    """Load gene pathogenicity scores for each patient."""
    df = pd.read_csv(scores_file, delimiter=',')
    
    # Get gene names from column headers
    genes = df.columns[1:].tolist()  # Skip patient_barcode column
    
    # Create a dictionary of patient scores
    patient_scores = {}
    for _, row in df.iterrows():
        patient_id = row['patient_barcode']
        gene_scores = {}
        for gene in genes:
            if not pd.isna(row[gene]):  # Only include non-NA scores
                gene_scores[gene] = float(row[gene])
        patient_scores[patient_id] = gene_scores
    
    return patient_scores

def aggregate_pathway_scores(pathways, patient_scores, aggregation_method='sum', score_threshold=0.5):
    """Calculate pathway scores for each patient.
    
    Args:
        pathways (dict): Dictionary mapping pathway names to sets of genes
        patient_scores (dict): Dictionary mapping patient IDs to gene scores
        aggregation_method (str): Method to aggregate gene scores ('sum', 'mean', 'median', 'max', or 'count_threshold')
        score_threshold (float): Minimum score threshold for including genes in calculations
    """
    pathway_scores = defaultdict(dict)
    
    for patient_id, scores in patient_scores.items():
        for pathway_name, pathway_genes in pathways.items():
            # Get scores for genes in this pathway that have scores
            pathway_gene_scores = [
                scores[gene] 
                for gene in pathway_genes 
                if gene in scores
            ]
            
            if pathway_gene_scores:  # Only calculate if we have scores
                if aggregation_method == 'sum':
                    # Only sum scores above the threshold
                    filtered_scores = [score for score in pathway_gene_scores if score >= score_threshold]
                    score = np.sum(filtered_scores) if filtered_scores else 0
                elif aggregation_method == 'count_threshold':
                    # Count number of scores above threshold
                    score = sum(1 for score in pathway_gene_scores if score >= score_threshold)
                elif aggregation_method == 'mean':
                    score = np.mean(pathway_gene_scores)
                elif aggregation_method == 'median':
                    score = np.median(pathway_gene_scores)
                elif aggregation_method == 'max':
                    score = np.max(pathway_gene_scores)
                else:
                    raise ValueError(f"Unknown aggregation method: {aggregation_method}")
                
                pathway_scores[pathway_name][patient_id] = score
    
    return pathway_scores

def create_pathway_score_matrix(pathway_scores):
    """Convert pathway scores dictionary to a DataFrame."""
    return pd.DataFrame(pathway_scores)

def main(pathogenicity_file, pathway_file, output_file, aggregation_method='sum', score_threshold=0.5):
    # Load data
    print("Loading pathway definitions...")
    pathways = load_pathway_genes(pathway_file)
    print(f"Loaded {len(pathways)} pathways")
    
    print("\nLoading pathogenicity scores...")
    patient_scores = load_pathogenicity_scores(pathogenicity_file)
    print(f"Loaded scores for {len(patient_scores)} patients")
    
    # Calculate pathway scores
    print("\nCalculating pathway scores...")
    print(f"Using {aggregation_method} method with score threshold of {score_threshold}")
    pathway_scores = aggregate_pathway_scores(
        pathways, 
        patient_scores, 
        aggregation_method,
        score_threshold
    )
    
    # Create and save results matrix
    results_df = create_pathway_score_matrix(pathway_scores)
    results_df.to_csv(output_file)
    print(f"\nSaved pathway scores to {output_file}")
    print(f"Final matrix shape: {results_df.shape}")
    
    # Print some basic statistics
    print("\nSummary statistics:")
    print(f"Number of pathways with scores: {results_df.shape[1]}")
    print(f"Number of patients with scores: {results_df.shape[0]}")
    
    return results_df

if __name__ == "__main__":
    # Example usage
    pathogenicity_file = "/Users/stanleychen/git/Melanoma/final_data/patient_maf_matrix_with_0s.csv"
    pathway_file = "/Users/stanleychen/git/Melanoma/data/gene_pathways/c2.cgp.v2024.1.Hs.symbols.gmt"
    output_file = "pathway_maf_scores_binary_above_cutoff.csv"
    
    results = main(
        pathogenicity_file,
        pathway_file,
        output_file,
        aggregation_method='count_threshold',  # Use new count_threshold method
        score_threshold=0.45
    )