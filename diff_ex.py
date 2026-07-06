import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

def load_and_prepare_data(expression_file, clinical_file):
    """
    Load and prepare the expression and clinical data
    """
    # Read the expression and clinical data
    expression_data = pd.read_csv(expression_file)
    clinical_data = pd.read_csv(clinical_file)
    
    # Set the first column (gene names) as index for expression data
    expression_data.set_index(expression_data.columns[0], inplace=True)
    
    # Create a mapping of patient to stage
    stage_map = dict(zip(clinical_data['patient_barcode'], 
                        clinical_data['tumor_stage']))
    
    return expression_data, stage_map

def perform_differential_expression(expression_data, stage_map):
    """
    Perform differential expression analysis between stage 1 and stage 2
    """
    results = []
    
    # For each gene
    for gene in expression_data.index:
        # Get expression values for stage 1 and stage 2
        stage1_expr = []
        stage2_expr = []
        
        # Group expression values by stage
        for patient in expression_data.columns:
            if patient in stage_map:
                expr_value = expression_data.loc[gene, patient]
                # Check for non-null and non-empty string values
                if pd.notnull(expr_value) and expr_value != '':
                    try:
                        value = float(expr_value)
                        if stage_map[patient] == 1.0:
                            stage1_expr.append(value)
                        elif stage_map[patient] == 2.0:
                            stage2_expr.append(value)
                    except ValueError:
                        continue
        
        # Only proceed if we have enough data for both stages
        if len(stage1_expr) >= 2 and len(stage2_expr) >= 2:  # Need at least 2 samples per group for t-test
            # Calculate statistics
            mean_stage1 = np.mean(stage1_expr)
            mean_stage2 = np.mean(stage2_expr)
            fold_change = abs(mean_stage2 - mean_stage1)
            
            # Check for variance in both groups
            var_stage1 = np.var(stage1_expr, ddof=1)
            var_stage2 = np.var(stage2_expr, ddof=1)
            
            # Only perform t-test if there's variance in the data
            if var_stage1 > 0 and var_stage2 > 0:
                try:
                    t_stat, p_value = stats.ttest_ind(stage1_expr, stage2_expr)
                    # Check if t_stat and p_value are finite
                    if np.isfinite(t_stat) and np.isfinite(p_value):
                        results.append({
                            'gene': gene,
                            'mean_stage1': mean_stage1,
                            'mean_stage2': mean_stage2,
                            'fold_change': fold_change,
                            'p_value': p_value,
                            't_statistic': t_stat,
                            'n_stage1': len(stage1_expr),
                            'n_stage2': len(stage2_expr),
                            'std_stage1': np.std(stage1_expr, ddof=1),
                            'std_stage2': np.std(stage2_expr, ddof=1)
                        })
                except:
                    # If t-test fails, still include the gene but with NA for p-value
                    results.append({
                        'gene': gene,
                        'mean_stage1': mean_stage1,
                        'mean_stage2': mean_stage2,
                        'fold_change': fold_change,
                        'p_value': np.nan,
                        't_statistic': np.nan,
                        'n_stage1': len(stage1_expr),
                        'n_stage2': len(stage2_expr),
                        'std_stage1': np.std(stage1_expr, ddof=1),
                        'std_stage2': np.std(stage2_expr, ddof=1)
                    })
    
    # Convert to DataFrame and sort by absolute fold change
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Sort by fold change, but only for rows with valid p-values
        results_df.sort_values('fold_change', ascending=False, inplace=True)
        
        # Calculate FDR only for non-NA p-values
        valid_pvals = results_df['p_value'].notna()
        if valid_pvals.any():
            _, q_values = fdrcorrection(results_df.loc[valid_pvals, 'p_value'])
            results_df.loc[valid_pvals, 'adj_p_value'] = q_values
    
    return results_df

def main():
    # File paths
    expression_file = '/Users/stanleychen/git/Melanoma/data/cleaned_RNA_expression_data.csv'
    clinical_file = '/Users/stanleychen/git/Melanoma/final_data/test_label.csv'
    
    # Load and prepare data
    print("Loading data...")
    expression_data, stage_map = load_and_prepare_data(expression_file, clinical_file)
    
    # Perform differential expression analysis
    print("Performing differential expression analysis...")
    results = perform_differential_expression(expression_data, stage_map)
    
    if len(results) > 0:
        # Get top 100 differentially expressed genes
        top_100 = results.head(100)
        
        # Save results
        top_100.to_csv('top_100_differential_expression_after_clean.csv', index=False)
        
        # Print summary
        print("\nTop 10 differentially expressed genes:")
        print(top_100.head(10)[['gene', 'fold_change', 'mean_stage1', 
                               'mean_stage2', 'p_value', 'adj_p_value']])
        print(f"\nTotal number of differentially expressed genes analyzed: {len(results)}")
    else:
        print("No significant results found. Please check your input data.")

if __name__ == "__main__":
    main()

