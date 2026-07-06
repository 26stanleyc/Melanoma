import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt
import os
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA

def filter_sparse_features(X, sparsity_threshold=0.95):
    """Remove features that are too sparse"""
    zero_ratio = (X == 0).mean()
    return zero_ratio < sparsity_threshold

def filter_collinear_features(X, threshold=0.9):
    """Remove features with high correlation."""
    corr_matrix = X.corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
    X_filtered_corr = X.drop(to_drop, axis=1)
    print(f"Dropped {len(to_drop)} highly correlated features.")
    return X_filtered_corr

def calculate_cindex(cph, X, T, E):
    """
    Calculate c-index correctly:
    - Higher hazard (risk) correlates with shorter survival time
    - Need to negate hazard so higher values align with longer survival for c-index
    """
    # Get feature columns only
    feature_cols = [col for col in X.columns if col not in ['survival_time_days', 'vital_status']]
    
    # Calculate hazard
    hazard = cph.predict_partial_hazard(X[feature_cols])
    
    # Negate hazard so higher values align with longer survival times
    return concordance_index(T, -hazard, E)

def plot_survival_curves(durations, events, predictions, title, output_path):
    """Create Kaplan-Meier survival plots stratified by binary risk groups (Low/High)"""
    # Create binary risk groups based on median predicted hazard
    risk_groups = pd.qcut(predictions, q=2, labels=['Low', 'High'])
    
    # Initialize the plot
    plt.figure(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    # Set colors similar to the image
    colors = {'Low': '#FF9999', 'High': '#66CCCC'}  # Light red and light blue
    
    # Plot each risk group with confidence intervals
    for group in ['Low', 'High']:
        mask = risk_groups == group
        n_patients = sum(mask)
        label = f"Risk = {group}"
        
        kmf.fit(
            durations[mask],
            events[mask],
            label=label
        )
        kmf.plot(
            ci_show=True,
            color=colors[group],
            alpha=0.7
        )
    
    plt.title(title)
    plt.xlabel('Survival Time (years)')
    plt.ylabel('Survival probability')
    plt.grid(True, alpha=0.3)
    
    # Add log-rank p-value
    from lifelines.statistics import logrank_test
    low_risk_mask = risk_groups == 'Low'
    high_risk_mask = risk_groups == 'High'
    p_value = logrank_test(
        durations[low_risk_mask], 
        durations[high_risk_mask],
        events[low_risk_mask],
        events[high_risk_mask]
    ).p_value
    plt.text(0.05, 0.05, f'P = {p_value:.4f}', 
             transform=plt.gca().transAxes)
    
    # Add number at risk table below the plot
    risk_counts = {}
    time_points = np.arange(0, int(max(durations)/365.25) + 1, 2)  # Every 2 years
    
    plt.grid(True, alpha=0.3)
    
    # Create a new axes for the table
    table_ax = plt.axes([0.1, 0.0, 0.8, 0.2])
    table_ax.axis('off')
    
    # Calculate patients at risk for each time point
    for group in ['Low', 'High']:
        mask = risk_groups == group
        risk_counts[group] = []
        for t in time_points:
            count = sum((durations[mask] >= t * 365.25))  # Convert years to days
            risk_counts[group].append(str(count))
    
    # Create the table
    rows = [risk_counts['Low'], risk_counts['High']]
    table = plt.table(
        cellText=rows,
        rowLabels=['-', '-'],
        colLabels=[f'{int(t)}' for t in time_points],
        loc='center',
        cellLoc='center'
    )
    
    # Add "Number of patients at risk" label
    plt.figtext(0.5, 0.15, 'Number of patients at risk', ha='center')
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def train_cox_model(features_file, label_file, output_dir, sparsity_threshold=0.95, 
                   correlation_threshold=0.9, penalizer=0.3, l1_ratio=0.5, cv=5):
    """Train Cox model with correct c-index calculation"""
    os.makedirs(output_dir, exist_ok=True)

    # Load and merge data
    print("Loading data...")
    features_df = pd.read_csv(features_file)
    labels_df = pd.read_csv(label_file)
    
    merged_df = pd.merge(features_df, 
                        labels_df[['patient_barcode', 'survival_time_days', 'vital_status']], 
                        on='patient_barcode', 
                        how='inner')
    
    print(f"Total samples after merging: {len(merged_df)}")

    # Prepare features and target
    X = merged_df.drop(['patient_barcode', 'survival_time_days', 'vital_status'], axis=1)
    merged_df['vital_status'] = (merged_df['vital_status'] == 'Dead').astype(int)
    T = merged_df['survival_time_days']
    E = merged_df['vital_status']

    # Filter features
    keep_features = filter_sparse_features(X, sparsity_threshold)
    X = X.loc[:, keep_features]
    print(f"Features after sparsity filtering: {X.shape}")
    
    X_filtered = filter_collinear_features(X, correlation_threshold)
    print(f"Features after correlation filtering: {X_filtered.shape}")

    # Cross-validation
    kf = KFold(n_splits=cv, shuffle=True, random_state=42)
    cv_c_indices = []

    for fold, (train_idx, val_idx) in enumerate(kf.split(X_filtered)):
        print(f"Fold {fold+1}/{cv}")
        
        # Split data
        X_train_fold = X_filtered.iloc[train_idx]
        X_val_fold = X_filtered.iloc[val_idx]
        T_train_fold = T.iloc[train_idx]
        T_val_fold = T.iloc[val_idx]
        E_train_fold = E.iloc[train_idx]
        E_val_fold = E.iloc[val_idx]

        # Scale features
        scaler = StandardScaler()
        X_train_fold_scaled = pd.DataFrame(
            scaler.fit_transform(X_train_fold),
            columns=X_train_fold.columns,
            index=X_train_fold.index
        )
        X_val_fold_scaled = pd.DataFrame(
            scaler.transform(X_val_fold),
            columns=X_val_fold.columns,
            index=X_val_fold.index
        )

        # Add survival data back for model fitting
        X_train_fold_scaled['survival_time_days'] = T_train_fold
        X_train_fold_scaled['vital_status'] = E_train_fold
        X_val_fold_scaled['survival_time_days'] = T_val_fold
        X_val_fold_scaled['vital_status'] = E_val_fold

        # Train model
        cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
        cph.fit(X_train_fold_scaled, 
                duration_col='survival_time_days', 
                event_col='vital_status')

        # Calculate c-index
        c_index_fold = calculate_cindex(cph, X_val_fold_scaled, T_val_fold, E_val_fold)
        cv_c_indices.append(c_index_fold)

    cv_c_indices = np.array(cv_c_indices)
    print(f"\nCross-validation C-index: {cv_c_indices.mean():.3f} (+/- {cv_c_indices.std():.3f})")

    # Final model training
    X_train, X_test, T_train, T_test, E_train, E_test = train_test_split(
        X_filtered, T, E, test_size=0.2, random_state=42)

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = pd.DataFrame(
        scaler.fit_transform(X_train),
        columns=X_train.columns,
        index=X_train.index
    )
    X_test_scaled = pd.DataFrame(
        scaler.transform(X_test),
        columns=X_test.columns,
        index=X_test.index
    )

    # Add survival data back
    X_train_scaled['survival_time_days'] = T_train
    X_train_scaled['vital_status'] = E_train
    X_test_scaled['survival_time_days'] = T_test
    X_test_scaled['vital_status'] = E_test

    # Train final model
    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
    cph.fit(X_train_scaled, duration_col='survival_time_days', event_col='vital_status')

    # Calculate final metrics
    train_cindex = calculate_cindex(cph, X_train_scaled, T_train, E_train)
    test_cindex = calculate_cindex(cph, X_test_scaled, T_test, E_test)

    # Generate predictions for survival plots
    feature_cols = [col for col in X_train_scaled.columns 
                   if col not in ['survival_time_days', 'vital_status']]
    train_predictions = cph.predict_partial_hazard(X_train_scaled[feature_cols])
    test_predictions = cph.predict_partial_hazard(X_test_scaled[feature_cols])

    # Create survival plots
    plot_survival_curves(
        T_train, E_train, train_predictions,
        'Training Set Survival Curves by Risk Group',
        os.path.join(output_dir, "train_survival_curves.png")
    )
    
    plot_survival_curves(
        T_test, E_test, test_predictions,
        'Test Set Survival Curves by Risk Group',
        os.path.join(output_dir, "test_survival_curves.png")
    )

    # Plot coefficients
    plt.figure(figsize=(12, 6))
    cph.plot()
    plt.title('Cox Model Coefficients')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cox_coefficients.png"))
    plt.close()

    print("\nFinal Results:")
    print(f"Training C-index: {train_cindex:.3f}")
    print(f"Test C-index: {test_cindex:.3f}")

    # Save results
    with open(os.path.join(output_dir, "cox_results.txt"), 'w') as f:
        f.write("Cox Proportional Hazards Model Results\n")
        f.write("=====================================\n\n")
        f.write(f"Data Statistics:\n")
        f.write(f"Original features: {len(keep_features)}\n")
        f.write(f"Features after filtering: {X_filtered.shape[1]}\n\n")
        f.write("Model Parameters:\n")
        f.write(f"Penalizer: {penalizer}\n")
        f.write(f"L1 ratio: {l1_ratio}\n\n")
        f.write("Metrics:\n")
        f.write(f"Training C-index: {train_cindex:.3f}\n")
        f.write(f"Test C-index: {test_cindex:.3f}\n")
        f.write(f"Cross-validation C-index: {cv_c_indices.mean():.3f} (+/- {cv_c_indices.std():.3f})\n\n")
        
        # Add top features by absolute coefficient value
        f.write("Top 10 Most Important Features:\n")
        coefficients = pd.Series(cph.params_)
        top_features = coefficients.abs().sort_values(ascending=False).head(10)
        for feature, coef in top_features.items():
            f.write(f"{feature}: {coefficients[feature]:.3f}\n")

    return cph, {
        'train_cindex': train_cindex,
        'test_cindex': test_cindex,
        'cv_cindex_mean': cv_c_indices.mean(),
        'cv_cindex_std': cv_c_indices.std(),
        'train_predictions': train_predictions,
        'test_predictions': test_predictions
    }

if __name__ == "__main__":
    # File paths
    features_file = "/Users/stanleychen/git/Melanoma/data/pathway_maf_scores_filtered_binary_above_cutoff.csv"
    label_file = "/Users/stanleychen/git/Melanoma/data/filtered_patient_data_all_v2.csv"
    output_dir = "/Users/stanleychen/git/Melanoma/BPSO/cox_results"
    
    # Train model
    cox_model, results = train_cox_model(
        features_file,
        label_file,
        output_dir,
        sparsity_threshold=0.9,
        correlation_threshold=0.65
    )
