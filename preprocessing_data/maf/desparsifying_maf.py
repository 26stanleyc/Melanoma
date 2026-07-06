import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from pathlib import Path
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

class CGPPathwayTransformer:
    def __init__(self, cgp_file_path: str = None):
        self.pathway_dict = self._load_cgp_pathways(cgp_file_path)
        
    def _load_cgp_pathways(self, file_path: str = None) -> Dict[str, List[str]]:
        """Load CGP gene sets from MSigDB GMT file"""
        pathway_dict = {}
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        pathway_name = parts[0]
                        genes = parts[2:]
                        pathway_dict[pathway_name] = genes
        except Exception as e:
            print(f"Error loading CGP pathways: {str(e)}")
            return {}
        return pathway_dict

    def _get_patient_from_sample(self, sample_id: str) -> str:
        """Extract patient ID from sample ID"""
        # Assuming format like TCGA-XX-YYYY-01A where TCGA-XX-YYYY is the patient
        parts = sample_id.split('-')
        if len(parts) >= 3:
            return '-'.join(parts[:3])
        return sample_id

    def filter_pathways(self, mutation_df: pd.DataFrame, 
                       min_mutations_per_pathway: int = 5,
                       min_pathway_frequency: float = 0.05) -> Dict[str, List[str]]:
        """
        Filter pathways based on mutation frequency and coverage
        
        Args:
            mutation_df: Mutation data
            min_mutations_per_pathway: Minimum mutations in pathway
            min_pathway_frequency: Minimum fraction of patients with pathway altered
        """
        filtered_pathways = {}
        unique_patients = mutation_df['patient_barcode'].apply(self._get_patient_from_sample).unique()
        total_patients = len(unique_patients)

        for pathway_name, genes in self.pathway_dict.items():
            # Get mutations in pathway genes
            pathway_mutations = mutation_df[mutation_df['gene'].isin(genes)]
            
            # Count unique patients with mutations in pathway
            affected_patients = pathway_mutations['patient_barcode'].apply(
                self._get_patient_from_sample).nunique()
            
            # Check if pathway meets criteria
            if (len(pathway_mutations) >= min_mutations_per_pathway and 
                affected_patients / total_patients >= min_pathway_frequency):
                filtered_pathways[pathway_name] = genes

        print(f"Filtered from {len(self.pathway_dict)} to {len(filtered_pathways)} pathways")
        return filtered_pathways

    def create_feature_matrix(self, mutation_df: pd.DataFrame,
                            min_mutations_per_pathway: int = 5,
                            min_pathway_frequency: float = 0.05,
                            variance_threshold: float = 0.01) -> Tuple[pd.DataFrame, List[str]]:
        """Create pathway-based feature matrix with dimensionality reduction"""
        
        # First filter pathways
        valid_pathways = self.filter_pathways(mutation_df, 
                                            min_mutations_per_pathway,
                                            min_pathway_frequency)
        
        # Get unique patients (not samples)
        patients = mutation_df['patient_barcode'].apply(self._get_patient_from_sample).unique()
        
        # Initialize feature matrices
        binary_matrix = pd.DataFrame(0, index=patients, columns=valid_pathways.keys())
        alpha_score_matrix = pd.DataFrame(0.0, index=patients, columns=valid_pathways.keys())
        
        # Process each patient
        for patient in patients:
            # Get all samples for this patient
            patient_samples = mutation_df[
                mutation_df['patient_barcode'].apply(self._get_patient_from_sample) == patient]
            
            # Aggregate mutations across all samples for this patient
            for pathway_name, genes in valid_pathways.items():
                pathway_mutations = patient_samples[patient_samples['gene'].isin(genes)]
                
                if len(pathway_mutations) > 0:
                    # Binary feature - pathway is altered
                    binary_matrix.loc[patient, pathway_name] = 1
                    
                    # Mean AlphaMissense score for pathway
                    valid_scores = pathway_mutations['AlphaMissense_Score'].dropna()
                    if len(valid_scores) > 0:
                        alpha_score_matrix.loc[patient, pathway_name] = valid_scores.mean()
        
        # Combine features
        feature_matrices = []
        
        # Add binary features with suffix
        binary_matrix.columns = [f"{col}_altered" for col in binary_matrix.columns]
        feature_matrices.append(binary_matrix)
        
        # Add AlphaMissense score features with suffix
        alpha_score_matrix.columns = [f"{col}_alpha_score" for col in alpha_score_matrix.columns]
        feature_matrices.append(alpha_score_matrix)
        
        # Combine all features
        combined_matrix = pd.concat(feature_matrices, axis=1)
        
        # Apply variance threshold to remove near-constant features
        selector = VarianceThreshold(threshold=variance_threshold)
        selected_features = selector.fit_transform(combined_matrix)
        
        # Create new dataframe with selected features
        selected_columns = combined_matrix.columns[selector.get_support()].tolist()
        final_matrix = pd.DataFrame(selected_features, 
                                  index=combined_matrix.index,
                                  columns=selected_columns)
        
        return final_matrix, selected_columns

    def add_global_features(self, feature_matrix: pd.DataFrame, 
                          mutation_df: pd.DataFrame) -> pd.DataFrame:
        """Add patient-level global features"""
        global_features = {}
        
        for patient in feature_matrix.index:
            # Get all samples for this patient
            patient_samples = mutation_df[
                mutation_df['patient_barcode'].apply(self._get_patient_from_sample) == patient]
            
            features = {
                'total_mutations': len(patient_samples),
                'mean_alpha_score': patient_samples['AlphaMissense_Score'].mean(),
                'total_pathogenic': len(patient_samples[
                    patient_samples['AlphaMissense_Class'] == 'likely_pathogenic'
                ])
            }
            
            global_features[patient] = features
            
        global_df = pd.DataFrame.from_dict(global_features, orient='index')
        enhanced_matrix = pd.concat([feature_matrix, global_df], axis=1)
        
        return enhanced_matrix

def prepare_dl_input(mutation_file: str, cgp_file: str, output_file: str = None,
                    min_mutations: int = 5, min_frequency: float = 0.05,
                    variance_threshold: float = 0.01) -> Tuple[np.ndarray, List[str]]:
    """Prepare mutation data for deep learning with dimensionality reduction"""
    
    # Read mutation data
    mutations = pd.read_csv(mutation_file)
    
    # Create transformer and generate features
    transformer = CGPPathwayTransformer(cgp_file)
    pathway_features, used_pathways = transformer.create_feature_matrix(
        mutations, min_mutations, min_frequency, variance_threshold
    )
    final_features = transformer.add_global_features(pathway_features, mutations)
    
    # Standardize features
    scaler = StandardScaler()
    scaled_features = pd.DataFrame(
        scaler.fit_transform(final_features),
        index=final_features.index,
        columns=final_features.columns
    )
    
    # Save if output file specified
    if output_file:
        scaled_features.to_csv(output_file)
        print(f"Saved feature matrix to {output_file}")
        
        # Save feature descriptions
        desc_file = output_file.rsplit('.', 1)[0] + '_features.txt'
        with open(desc_file, 'w') as f:
            f.write("Feature Matrix Description\n")
            f.write(f"Number of patients: {scaled_features.shape[0]}\n")
            f.write(f"Number of features: {scaled_features.shape[1]}\n\n")
            f.write("Features:\n")
            for i, feature in enumerate(scaled_features.columns):
                f.write(f"{i+1}. {feature}\n")
        print(f"Saved feature descriptions to {desc_file}")
    
    return scaled_features.values, scaled_features.columns.tolist()

if __name__ == "__main__":
    mutation_file = "/Users/stanleychen/git/Melanoma/final_data/mutations_with_alphamissense.csv"
    cgp_file = "/Users/stanleychen/git/Melanoma/data/gene_pathways/c2.cgp.v2024.1.Hs.symbols.gmt"
    output_file = "maf_features.csv"
    
    X, features = prepare_dl_input(
        mutation_file, 
        cgp_file, 
        output_file,
        min_mutations=5,          # Minimum mutations in pathway
        min_frequency=0.05,       # Pathway must be altered in 5% of patients
        variance_threshold=0.05   # Remove features with variance < 0.01
    )
    
    print(f"Feature matrix shape: {X.shape}")
    print(f"Number of features: {len(features)}")