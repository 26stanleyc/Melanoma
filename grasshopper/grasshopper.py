import numpy as np
from sklearn.preprocessing import StandardScaler
from typing import List, Tuple, Dict
import pandas as pd
import json
import os
from datetime import datetime
from scipy import stats

# Configure paths
TUMOR_FILE = "/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/filtered_rna_seq-R.csv"
CONTROL_FILE = "/Users/stanleychen/git/Melanoma/final_data/control_folder_v2/filtered_expression_log2_tpm-R.csv"
OUTPUT_DIR = "/Users/stanleychen/git/Melanoma/grasshopper/results"


class GOL2_2T:
    def __init__(
        self,
        population_size: int = 50,
        max_iterations: int = 100,
        c_init: float = 2.0,
        g: float = 0.1,
        U: float = 1.0,
        l2_alpha: float = 0.01
    ):
        self.population_size = population_size
        self.max_iterations = max_iterations
        self.c_init = c_init
        self.g = g
        self.U = U
        self.l2_alpha = l2_alpha
        self.best_solution = None
        self.best_fitness = float('-inf')
        self.fitness_history = []

    def _calculate_l2_regularization(self, x: np.ndarray) -> float:
        """Calculate L2 regularization term."""
        return np.sum(x ** 2)
        
    def load_data(
        self,
        tumor_filepath: str,
        control_filepath: str
    ) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """Load tumor and control RNA-seq data."""
        tumor_df = pd.read_csv(tumor_filepath, index_col=0)
        control_df = pd.read_csv(control_filepath, index_col=0)
        
        common_genes = sorted(set(tumor_df.index).intersection(set(control_df.index)))
        
        if len(common_genes) == 0:
            raise ValueError("No common genes found between datasets")
            
        print(f"Found {len(common_genes)} common genes")
        print(f"Tumor samples: {tumor_df.shape[1]}")
        print(f"Control samples: {control_df.shape[1]}")
        
        tumor_data = tumor_df.loc[common_genes].values.astype(float)
        control_data = control_df.loc[common_genes].values.astype(float)
        
        return tumor_data.T, control_data.T, common_genes
    
    def _calculate_social_interaction(
        self,
        x_i: np.ndarray,
        population: np.ndarray,
        c: float
    ) -> np.ndarray:
        """Calculate social interaction component."""
        N = population.shape[0]
        diff_sum = np.sum(population, axis=0) - N * x_i
        return c * (diff_sum / N)
    
    def _calculate_gravity(self, x_i: np.ndarray) -> np.ndarray:
        """Calculate gravity component."""
        o = 0.5
        return -self.g * (x_i - o)
    
    def _calculate_wind(
        self,
        x_i: np.ndarray,
        x_best: np.ndarray
    ) -> np.ndarray:
        """Calculate wind component."""
        return self.U * (x_best - x_i)

    def _calculate_differential_metrics(
        self,
        tumor_data: np.ndarray,
        control_data: np.ndarray,
        selected: np.ndarray
    ) -> Dict[str, float]:
        """Calculate differential expression metrics safely."""
        if not np.any(selected):
            return {
                'fold_change': 0,
                'ttest_score': 0,
                'effect_size': 0,
                'variance_ratio': 0
            }
        
        # Get selected gene data
        tumor_selected = tumor_data[:, selected]
        control_selected = control_data[:, selected]
        
        # Calculate mean expression
        tumor_means = np.maximum(np.mean(tumor_selected, axis=0), 1e-10)  # Avoid log(0)
        control_means = np.maximum(np.mean(control_selected, axis=0), 1e-10)
        
        # Safe log2 fold change calculation
        fold_changes = np.log2((tumor_means + 1) / (control_means + 1))
        
        # T-test with safety checks
        ttest_scores = []
        for i in range(tumor_selected.shape[1]):
            try:
                t_stat, _ = stats.ttest_ind(
                    tumor_selected[:, i],
                    control_selected[:, i],
                    equal_var=False  # Welch's t-test
                )
                ttest_scores.append(abs(t_stat))
            except:
                ttest_scores.append(0)
        
        # Effect size calculation with safety checks
        effect_sizes = []
        for i in range(tumor_selected.shape[1]):
            tumor_var = np.var(tumor_selected[:, i])
            control_var = np.var(control_selected[:, i])
            pooled_std = np.sqrt((tumor_var + control_var) / 2 + 1e-10)
            effect_size = abs(
                (np.mean(tumor_selected[:, i]) - np.mean(control_selected[:, i])) / 
                pooled_std
            )
            effect_sizes.append(effect_size)
        
        # Safe variance ratio calculation
        tumor_vars = np.var(tumor_selected, axis=0) + 1e-10
        control_vars = np.var(control_selected, axis=0) + 1e-10
        variance_ratios = np.log2(tumor_vars / control_vars)
        
        return {
            'fold_change': float(np.mean(np.abs(fold_changes))),
            'ttest_score': float(np.mean(ttest_scores)),
            'effect_size': float(np.mean(effect_sizes)),
            'variance_ratio': float(np.mean(np.abs(variance_ratios)))
        }
    
    def _evaluate_fitness(
        self,
        position: np.ndarray,
        tumor_data: np.ndarray,
        control_data: np.ndarray
    ) -> float:
        """Evaluate fitness based on tumor vs control differences."""
        selected = position > 0.5
        if not np.any(selected):
            return float('-inf')
        
        metrics = self._calculate_differential_metrics(
            tumor_data,
            control_data,
            selected
        )
        
        # Combine metrics with weights
        fitness = (
            5.0 * metrics['fold_change'] +
            2.0 * metrics['effect_size'] +
            1.0 * metrics['ttest_score'] +
            1.0 * metrics['variance_ratio'] -
            self.l2_alpha * self._calculate_l2_regularization(position)
        )
        
        return float(fitness)
    
    def fit(
        self,
        tumor_data: np.ndarray,
        control_data: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """Fit the model using tumor and control data."""
        try:
            n_features = tumor_data.shape[1]
            population = np.random.uniform(0, 1, (self.population_size, n_features))
            
            self.best_solution = population[0].copy()
            self.best_fitness = self._evaluate_fitness(
                self.best_solution,
                tumor_data,
                control_data
            )
            self.fitness_history = [self.best_fitness]
            
            for iteration in range(self.max_iterations):
                c = self.c_init * (1 - iteration / self.max_iterations)
                
                for i in range(self.population_size):
                    current_position = population[i].copy()
                    
                    S = self._calculate_social_interaction(
                        current_position,
                        population,
                        c
                    )
                    G = self._calculate_gravity(current_position)
                    A = self._calculate_wind(current_position, self.best_solution)
                    
                    new_position = current_position + S + G + A
                    population[i] = np.clip(new_position, 0, 1)
                    
                    fitness = self._evaluate_fitness(
                        population[i],
                        tumor_data,
                        control_data
                    )
                    
                    if fitness > self.best_fitness:
                        self.best_solution = population[i].copy()
                        self.best_fitness = fitness
                
                if (iteration + 1) % 10 == 0:
                    n_selected = np.sum(self.best_solution > 0.5)
                    print(f"Iteration {iteration + 1}/{self.max_iterations}, "
                          f"Best Fitness: {self.best_fitness:.4f}, "
                          f"Selected Features: {n_selected}")
            
            feature_mask = self.best_solution > 0.5
            
            return {
                'best_solution': self.best_solution,
                'feature_mask': feature_mask,
                'best_fitness': self.best_fitness,
                'fitness_history': self.fitness_history
            }
            
        except Exception as e:
            print(f"Error in fit method: {str(e)}")
            raise

def main():
    try:
        model = GOL2_2T()
        
        print("Loading data...")
        tumor_data, control_data, gene_names = model.load_data(
            TUMOR_FILE,
            CONTROL_FILE
        )
        
        print("\nPreprocessing data...")
        scaler = StandardScaler()
        all_samples = np.vstack([tumor_data, control_data])
        scaler.fit(all_samples)
        
        tumor_scaled = scaler.transform(tumor_data)
        control_scaled = scaler.transform(control_data)
        
        print(f"\nData shapes:")
        print(f"Tumor data: {tumor_scaled.shape}")
        print(f"Control data: {control_scaled.shape}")
        
        print("\nStarting feature selection...")
        results = model.fit(tumor_scaled, control_scaled)
        
        # Save results
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        selected_indices = np.where(results['feature_mask'])[0]
        selected_genes = [gene_names[i] for i in selected_indices]
        
        save_dict = {
            "timestamp": timestamp,
            "selected_genes": selected_genes,
            "num_selected_genes": len(selected_genes),
            "best_fitness": float(results['best_fitness'])
        }
        
        output_file = os.path.join(OUTPUT_DIR, f"differential_genes_{timestamp}.json")
        with open(output_file, 'w') as f:
            json.dump(save_dict, f, indent=4)
        
        print(f"\nSelected {len(selected_genes)} genes")
        print(f"Results saved to: {output_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()