import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index
from tensorflow.keras import layers, models, callbacks, regularizers
import matplotlib.pyplot as plt
import os

def filter_sparse_features(X, sparsity_threshold=0.95):
    """Remove features that are too sparse"""
    zero_ratio = (X == 0).mean()
    return zero_ratio < sparsity_threshold

def create_model(input_dim):
    """Create neural network model for sparse data"""
    model = models.Sequential([
        # First layer with strong L1 regularization
        layers.Dense(64, activation='relu', input_dim=input_dim,
                    kernel_regularizer=regularizers.l1_l2(l1=0.01, l2=0.01)),
        layers.BatchNormalization(),
        layers.Dropout(0.4),
        
        # Hidden layers with reducing dimensions
        layers.Dense(32, activation='relu',
                    kernel_regularizer=regularizers.l1_l2(l1=0.01, l2=0.01)),
        layers.BatchNormalization(),
        layers.Dropout(0.3),
        
        # Output layer
        layers.Dense(1, activation='linear')
    ])
    
    model.compile(
        optimizer='adam',
        loss='mse',
        metrics=['mae']
    )
    return model

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
    table_ax = plt.axes([0.1, 0.0, 0.8, 0.2])  # Adjust position as needed
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
        rowLabels=['-', '-'],  # Matching the image style
        colLabels=[f'{int(t)}' for t in time_points],
        loc='center',
        cellLoc='center'
    )
    
    # Add "Number of patients at risk" label
    plt.figtext(0.5, 0.15, 'Number of patients at risk', ha='center')
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def train_survival_nn(features_file, label_file, output_dir, sparsity_threshold=0.95):
    """Train Neural Network for sparse survival data"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("Loading data...")
    features_df = pd.read_csv(features_file)
    labels_df = pd.read_csv(label_file)
    
    # Merge datasets
    merged_df = pd.merge(features_df, 
                        labels_df[['patient_barcode', 'survival_time_days', 'vital_status']], 
                        on='patient_barcode',
                        how='inner')
    
    print(f"Total samples after merging: {len(merged_df)}")
    
    # Prepare features and target
    X = merged_df.drop(['patient_barcode', 'survival_time_days', 'vital_status'], axis=1)
    y = merged_df['survival_time_days']
    events = (merged_df['vital_status'] == 'Dead').astype(int)
    
    # Analyze sparsity
    sparsity = (X == 0).mean().mean()
    print(f"\nOverall data sparsity: {sparsity:.3f}")
    
    # Filter sparse features
    keep_features = filter_sparse_features(X, sparsity_threshold)
    X = X.loc[:, keep_features]
    print(f"Features retained after sparsity filtering: {X.shape[1]}")
    
    # Split data
    X_train, X_test, y_train, y_test, events_train, events_test = train_test_split(
        X, y, events, test_size=0.2, random_state=42
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Create and train model
    print("\nTraining neural network...")
    model = create_model(input_dim=X_train_scaled.shape[1])
    
    # Callbacks
    early_stopping = callbacks.EarlyStopping(
        monitor='val_loss',
        patience=30,
        restore_best_weights=True
    )
    
    reduce_lr = callbacks.ReduceLROnPlateau(
        monitor='val_loss',
        factor=0.5,
        patience=10,
        min_lr=0.0001
    )
    
    # Train with balanced class weights
    history = model.fit(
        X_train_scaled, y_train,
        validation_split=0.2,
        epochs=200,
        batch_size=32,
        callbacks=[early_stopping, reduce_lr],
        verbose=1
    )
    
    # Predictions
    y_pred_train = model.predict(X_train_scaled).flatten()
    y_pred_test = model.predict(X_test_scaled).flatten()
    
    # Calculate metrics
    train_mse = mean_squared_error(y_train, y_pred_train)
    test_mse = mean_squared_error(y_test, y_pred_test)
    train_r2 = r2_score(y_train, y_pred_train)
    test_r2 = r2_score(y_test, y_pred_test)
    train_cindex = concordance_index(y_train, y_pred_train, events_train)
    test_cindex = concordance_index(y_test, y_pred_test, events_test)
    
    print("\nFinal Results:")
    print(f"Training MSE: {train_mse:.2f}")
    print(f"Test MSE: {test_mse:.2f}")
    print(f"Training R²: {train_r2:.3f}")
    print(f"Test R²: {test_r2:.3f}")
    print(f"Training C-index: {train_cindex:.3f}")
    print(f"Test C-index: {test_cindex:.3f}")
    
    # Plot training history
    plt.figure(figsize=(10, 5))
    plt.plot(history.history['loss'], label='Training Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training History')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "nn_training_history.png"))
    plt.close()
    
    # Create survival plots
    plot_survival_curves(
        y_train, events_train, y_pred_train,
        'Training Set Survival Curves by Risk Group',
        os.path.join(output_dir, "train_survival_curves.png")
    )
    
    plot_survival_curves(
        y_test, events_test, y_pred_test,
        'Test Set Survival Curves by Risk Group',
        os.path.join(output_dir, "test_survival_curves.png")
    )
    
    # Save results
    with open(os.path.join(output_dir, "nn_results.txt"), 'w') as f:
        f.write("Neural Network Results for Sparse Data\n")
        f.write("====================================\n\n")
        f.write(f"Data Statistics:\n")
        f.write(f"Overall sparsity: {sparsity:.3f}\n")
        f.write(f"Original features: {len(keep_features)}\n")
        f.write(f"Features after filtering: {X.shape[1]}\n\n")
        f.write("Metrics:\n")
        f.write(f"Training MSE: {train_mse:.2f}\n")
        f.write(f"Test MSE: {test_mse:.2f}\n")
        f.write(f"Training R²: {train_r2:.3f}\n")
        f.write(f"Test R²: {test_r2:.3f}\n")
        f.write(f"Training C-index: {train_cindex:.3f}\n")
        f.write(f"Test C-index: {test_cindex:.3f}\n")
    
    return model, history, {
        'train_mse': train_mse,
        'test_mse': test_mse,
        'train_r2': train_r2,
        'test_r2': test_r2,
        'train_cindex': train_cindex,
        'test_cindex': test_cindex
    }

if __name__ == "__main__":
    # File paths
    features_file = "/Users/stanleychen/git/Melanoma/data/pathway_maf_scores_filtered_binary_above_cutoff.csv"
    label_file = "/Users/stanleychen/git/Melanoma/data/patient_survival_analysis_filtered_keep_negs.csv"
    output_dir = "/Users/stanleychen/git/Melanoma/BPSO/nn_results"
    
    # Train model
    model, history, results = train_survival_nn(
        features_file, 
        label_file, 
        output_dir,
        sparsity_threshold=0.95
    )