import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import os
from datetime import datetime

class BPSO_DT:
    def __init__(self, n_particles, n_iterations, target_features=500, feature_tolerance=50,
                 alpha=0.5, w=0.9, c1=2, c2=2):
        self.n_particles = n_particles
        self.n_iterations = n_iterations
        self.target_features = target_features
        self.feature_tolerance = feature_tolerance  # Allowable deviation from target
        self.alpha = alpha
        self.w = w
        self.c1 = c1
        self.c2 = c2
        self.history = []
        
    def _sigmoid(self, x):
        return 1 / (1 + np.exp(-x))
    
    def _calculate_fitness(self, X, y, position):
        """Calculate fitness with hard constraints on feature count"""
        n_selected = np.sum(position)
        
        # Return 0 fitness if feature count is outside acceptable range
        min_features = max(1, self.target_features - self.feature_tolerance)
        max_features = self.target_features + self.feature_tolerance
        
        if n_selected < min_features or n_selected > max_features:
            return 0.0
            
        selected_features = X[:, position == 1]
        
        X_train, X_test, y_train, y_test = train_test_split(
            selected_features, y, test_size=0.2, random_state=42
        )
        
        dt = DecisionTreeRegressor(random_state=42)
        dt.fit(X_train, y_train)
        
        y_pred = dt.predict(X_test)
        r2 = r2_score(y_test, y_pred)
        
        # Normalize feature count penalty between [0, 1]
        feature_diff = abs(n_selected - self.target_features)
        feature_penalty = 1 - (feature_diff / self.feature_tolerance)
        
        # Final fitness combines R² and feature count penalty
        fitness = self.alpha * r2 + (1 - self.alpha) * feature_penalty
        
        return fitness
    
    def _initialize_population(self, n_features):
        """Initialize population with controlled feature counts"""
        positions = np.zeros((self.n_particles, n_features))
        
        for i in range(self.n_particles):
            # Randomly select number of features within tolerance range
            n_selected = np.random.randint(
                self.target_features - self.feature_tolerance,
                self.target_features + self.feature_tolerance + 1
            )
            
            # Randomly select indices
            selected_indices = np.random.choice(
                n_features, size=n_selected, replace=False
            )
            positions[i, selected_indices] = 1
            
        return positions
    
    def fit(self, X, y, feature_names=None):
        n_features = X.shape[1]
        if feature_names is None:
            feature_names = [f'Feature_{i}' for i in range(n_features)]
        
        self.feature_names_ = feature_names
        
        # Initialize with controlled feature counts
        positions = self._initialize_population(n_features)
        velocities = np.random.randn(self.n_particles, n_features) * 0.1  # Reduced initial velocity
        
        pbest = positions.copy()
        pbest_fitness = np.array([self._calculate_fitness(X, y, pos) for pos in positions])
        
        gbest_idx = np.argmax(pbest_fitness)
        gbest = pbest[gbest_idx].copy()
        gbest_fitness = pbest_fitness[gbest_idx]
        
        # Adaptive inertia weight
        w_max = 0.9
        w_min = 0.4
        
        for iteration in range(self.n_iterations):
            # Update inertia weight
            w = w_max - (w_max - w_min) * iteration / self.n_iterations
            
            r1, r2 = np.random.rand(2)
            velocities = (w * velocities + 
                        self.c1 * r1 * (pbest - positions) +
                        self.c2 * r2 * (gbest - positions))
            
            # Clip velocities to prevent extreme values
            velocities = np.clip(velocities, -4, 4)
            
            prob = self._sigmoid(velocities)
            
            # Controlled update of positions
            new_positions = (prob > np.random.rand(self.n_particles, n_features)).astype(int)
            
            # Force feature count within bounds for some particles
            for i in range(self.n_particles):
                n_selected = np.sum(new_positions[i])
                if n_selected < self.target_features - self.feature_tolerance or \
                   n_selected > self.target_features + self.feature_tolerance:
                    if np.random.rand() < 0.5:  # 50% chance to force valid solution
                        n_desired = np.random.randint(
                            self.target_features - self.feature_tolerance,
                            self.target_features + self.feature_tolerance + 1
                        )
                        if n_selected > n_desired:
                            # Randomly remove features
                            ones = np.where(new_positions[i] == 1)[0]
                            to_remove = np.random.choice(ones, size=n_selected - n_desired, replace=False)
                            new_positions[i, to_remove] = 0
                        else:
                            # Randomly add features
                            zeros = np.where(new_positions[i] == 0)[0]
                            to_add = np.random.choice(zeros, size=n_desired - n_selected, replace=False)
                            new_positions[i, to_add] = 1
            
            positions = new_positions
            fitness = np.array([self._calculate_fitness(X, y, pos) for pos in positions])
            
            improved = fitness > pbest_fitness
            pbest[improved] = positions[improved]
            pbest_fitness[improved] = fitness[improved]
            
            best_idx = np.argmax(pbest_fitness)
            if pbest_fitness[best_idx] > gbest_fitness:
                gbest = pbest[best_idx].copy()
                gbest_fitness = pbest_fitness[best_idx]
            
            self.history.append({
                'iteration': iteration + 1,
                'best_fitness': gbest_fitness,
                'n_selected_features': np.sum(gbest)
            })
            
            print(f"Iteration {iteration + 1}/{self.n_iterations}, "
                  f"Best fitness: {gbest_fitness:.4f}, "
                  f"Selected features: {np.sum(gbest)}/{n_features} "
                  f"(target: {self.target_features}±{self.feature_tolerance})")
        
        self.best_position_ = gbest
        self.best_fitness_ = gbest_fitness
        
        return self
    
    def transform(self, X):
        return X[:, self.best_position_ == 1]
    
    def fit_transform(self, X, y, feature_names=None):
        self.fit(X, y, feature_names)
        return self.transform(X)
    
    def save_results(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # Save selected features with gene names
        selected_features_df = pd.DataFrame({
            'Gene': np.array(self.feature_names_)[self.best_position_ == 1],
            'Selected': True
        })
        selected_features_df.to_csv(
            os.path.join(output_dir, f'selected_genes_{timestamp}.csv'),
            index=False
        )
        
        # Save optimization history
        history_df = pd.DataFrame(self.history)
        history_df.to_csv(
            os.path.join(output_dir, f'optimization_history_{timestamp}.csv'),
            index=False
        )
        
        # Save summary with gene-specific terminology
        with open(os.path.join(output_dir, f'results_summary_{timestamp}.txt'), 'w') as f:
            f.write("BPSO-DT Gene Selection Results\n")
            f.write("============================\n\n")
            f.write(f"Parameters:\n")
            f.write(f"- Number of particles: {self.n_particles}\n")
            f.write(f"- Number of iterations: {self.n_iterations}\n")
            f.write(f"- Target number of genes: {self.target_features}\n")
            f.write(f"- Feature tolerance: ±{self.feature_tolerance}\n")
            f.write(f"- Alpha: {self.alpha}\n")
            f.write(f"- Inertia weight (w): {self.w}\n")
            f.write(f"- Cognitive coefficient (c1): {self.c1}\n")
            f.write(f"- Social coefficient (c2): {self.c2}\n\n")
            f.write(f"Results:\n")
            f.write(f"- Best fitness achieved: {self.best_fitness_:.4f}\n")
            f.write(f"- Number of selected genes: {np.sum(self.best_position_)}\n")
            f.write(f"- Total number of genes: {len(self.best_position_)}\n")
            f.write(f"- Deviation from target: {abs(np.sum(self.best_position_) - self.target_features)}\n")

def main():
    # File paths
    expression_file = '/Users/stanleychen/git/Melanoma/final_data/cleaned_RNA-seq_data.csv'  # Update with your file path
    label_file = '/Users/stanleychen/git/Melanoma/final_data/filtered_patient_data_all_v2.csv'  # Update with your file path
    output_dir = '/Users/stanleychen/git/Melanoma/BPSO/results/feature_selection'  # Update with your output directory

    # Load gene expression data
    print("Loading gene expression data...")
    expression_df = pd.read_csv(expression_file, index_col=0)
    
    # Load patient survival data
    print("Loading patient data...")
    label_df = pd.read_csv(label_file)
    
    print(f"Total genes: {len(expression_df.index)}")
    print(f"Total patients in expression data: {len(expression_df.columns)}")
    print(f"Total patients in label data: {len(label_df)}")
    
    # Extract patient IDs from expression data columns
    patient_ids = expression_df.columns
    
    # Match patients with survival data
    matched_patients = label_df[label_df['patient_barcode'].isin(patient_ids)]
    print(f"Matched patients: {len(matched_patients)}")
    
    # Prepare X matrix (genes as features) and y vector (survival times)
    X_matrix = expression_df[matched_patients['patient_barcode']].values.T
    y_vector = matched_patients['survival_time_days'].values
    gene_names = expression_df.index.tolist()
    
    # Initialize and run BPSO-DT
    print("\nInitializing BPSO-DT...")
    bpso = BPSO_DT(
        n_particles=30,
        n_iterations=100,
        target_features=500,
        feature_tolerance=50,  # Allow ±50 genes from target
        alpha=0.5
    )
    
    # Run feature selection
    print("\nRunning gene selection...")
    bpso.fit(X_matrix, y_vector, gene_names)
    
    # Save results
    print("\nSaving results...")
    bpso.save_results(output_dir)
    
    print("\nGene selection completed!")
    print(f"Results saved to: {output_dir}")
    print(f"Number of selected genes: {np.sum(bpso.best_position_)}")
    print(f"Best fitness achieved: {bpso.best_fitness_:.4f}")

if __name__ == "__main__":
    main()

