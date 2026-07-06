import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from torch.utils.data import Dataset, DataLoader
from lifelines.utils import concordance_index
from lifelines import KaplanMeierFitter
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, KFold
from sklearn.calibration import calibration_curve
from sksurv.metrics import cumulative_dynamic_auc
import warnings
warnings.filterwarnings('ignore')

def filter_sparse_features(X, sparsity_threshold=0.95):
    """Remove features that are too sparse"""
    zero_ratio = (X == 0).mean()
    return X.loc[:, zero_ratio < sparsity_threshold]

def process_and_align_data(mutation_file, mirna_file, embedding_file, id_mapping_file, labels_file):
    """
    Align data across all modalities using patient IDs with sparsity reduction
    """
    print("Loading and processing data...")
    
    # Load ID mappings
    id_mappings = pd.read_csv(id_mapping_file, sep='\t', names=['uuid', 'patient_id'])
    id_mappings = id_mappings.drop_duplicates(subset='patient_id', keep='first')
    print(f"Loaded {len(id_mappings)} unique ID mappings")
    
    # Load labels
    labels = pd.read_csv(labels_file)
    labels = labels.drop_duplicates(subset='patient_barcode', keep='first')
    valid_patients = set(labels['patient_barcode'])
    print(f"Loaded {len(valid_patients)} unique patients from labels")
    
    # Load and process mutation data with sparsity reduction
    mutations = pd.read_csv(mutation_file)
    mutations = mutations.drop_duplicates(subset='patient_barcode', keep='first')
    mutations = mutations[mutations['patient_barcode'].isin(valid_patients)]
    mutation_features = mutations.set_index('patient_barcode')
    
    # Apply sparsity filtering to mutation features
    print("Original mutation features:", mutation_features.shape[1])
    mutation_features = filter_sparse_features(mutation_features, sparsity_threshold=0.95)
    print("Mutation features after sparsity filtering:", mutation_features.shape[1])
    
    # Load miRNA data
    mirna = pd.read_csv(mirna_file, index_col=0)
    mirna=mirna.T
    mirna = mirna[~mirna.index.duplicated(keep='first')]
    mirna = mirna[mirna.index.isin(valid_patients)]
    print(f"Loaded miRNA data: {mirna.shape}")
    
    # Load embeddings
    embeddings = np.load(embedding_file)
    print(f"Loaded embeddings with shape: {embeddings.shape}")
    
    # Match embeddings with patient IDs
    matched_indices = []
    matched_patients = []
    
    for idx, (_, row) in enumerate(id_mappings.iterrows()):
        patient_id = row['patient_id']
        if patient_id in valid_patients and patient_id not in set(matched_patients):
            matched_indices.append(idx)
            matched_patients.append(patient_id)
    
    if not matched_indices:
        raise ValueError("No matching patients found between embeddings and labels!")
    
    embeddings = embeddings[matched_indices]
    embeddings_df = pd.DataFrame(
        embeddings,
        index=matched_patients,
        columns=[f'emb_{i}' for i in range(embeddings.shape[1])]
    )
    
    # Find common patients across all modalities
    common_patients = set(mutation_features.index) & set(mirna.index) & set(embeddings_df.index)
    print(f"Found {len(common_patients)} common patients across all modalities")
    
    if len(common_patients) == 0:
        raise ValueError("No common patients found across all modalities!")
    
    # Filter all data to common patients
    mutation_features = mutation_features[mutation_features.index.isin(common_patients)]
    mirna = mirna[mirna.index.isin(common_patients)]
    embeddings_df = embeddings_df[embeddings_df.index.isin(common_patients)]
    labels = labels[labels['patient_barcode'].isin(common_patients)]
    
    # Convert survival times and events
    labels['event'] = (labels['vital_status'] == 'Dead').astype(int)
    
    # Scale numerical features
    scaler = StandardScaler()
    mirna_scaled = pd.DataFrame(
        scaler.fit_transform(mirna),
        index=mirna.index,
        columns=mirna.columns
    )
    embeddings_scaled = pd.DataFrame(
        scaler.fit_transform(embeddings_df),
        index=embeddings_df.index,
        columns=embeddings_df.columns
    )
    mutation_scaled = pd.DataFrame(
        scaler.fit_transform(mutation_features),
        index=mutation_features.index,
        columns=mutation_features.columns
    )
    
    # Sort all dataframes by patient ID for consistency
    common_patients = sorted(list(common_patients))
    mutation_scaled = mutation_scaled.reindex(index=common_patients)
    mirna_scaled = mirna_scaled.reindex(index=common_patients)
    embeddings_scaled = embeddings_scaled.reindex(index=common_patients)
    labels = labels.set_index('patient_barcode').reindex(index=common_patients)
    
    return mutation_scaled, mirna_scaled, embeddings_scaled, labels

class MultiModalSurvivalDataset(Dataset):
    def __init__(self, mutation_data, mirna_data, embedding_data, labels):
        self.mutation_data = torch.FloatTensor(mutation_data.values)
        self.mirna_data = torch.FloatTensor(mirna_data.values)
        self.embedding_data = torch.FloatTensor(embedding_data.values)
        self.times = torch.FloatTensor(labels['survival_time_days'].values)
        self.events = torch.FloatTensor(labels['event'].values)
        
    def __len__(self):
        return len(self.times)
    
    def __getitem__(self, idx):
        return (self.mutation_data[idx], 
                self.mirna_data[idx], 
                self.embedding_data[idx]), \
               self.times[idx], \
               self.events[idx]

class MultiModalDeepSurv(nn.Module):
    def __init__(self, input_dims, hidden_dim=64, dropout=0.4):
        super(MultiModalDeepSurv, self).__init__()
        
        self.mutation_encoder = nn.Sequential(
            nn.Linear(input_dims[0], hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.BatchNorm1d(hidden_dim//2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        self.mirna_encoder = nn.Sequential(
            nn.Linear(input_dims[1], hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.BatchNorm1d(hidden_dim//2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        self.embedding_encoder = nn.Sequential(
            nn.Linear(input_dims[2], hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.BatchNorm1d(hidden_dim//2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        total_fusion_dim = (hidden_dim//2) * 3
        self.fusion_net = nn.Sequential(
            nn.Linear(total_fusion_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim//2),
            nn.BatchNorm1d(hidden_dim//2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim//2, 1)
        )
        
        self._init_weights()
    
    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                nn.init.zeros_(m.bias)
    
    def forward(self, x):
        mutation_data, mirna_data, embedding_data = x
        
        mutation_features = self.mutation_encoder(mutation_data)
        mirna_features = self.mirna_encoder(mirna_data)
        embedding_features = self.embedding_encoder(embedding_data)
        
        fused = torch.cat([mutation_features, mirna_features, embedding_features], dim=1)
        return self.fusion_net(fused)

def cox_loss(risk_pred, y_time, y_event):
    order = torch.argsort(y_time, descending=True)
    risk_pred = risk_pred[order]
    y_event = y_event[order]
    
    log_risk = torch.log(torch.cumsum(torch.exp(risk_pred), dim=0))
    likelihood = risk_pred - log_risk
    uncensored_likelihood = likelihood * y_event
    
    return -torch.mean(uncensored_likelihood)

def create_nomogram(model, feature_names, device, scaler=None):
    """Create a nomogram visualization for the model"""
    plt.figure(figsize=(12, 8))
    
    with torch.no_grad():
        mutation_weights = torch.norm(model.mutation_encoder[0].weight.data, dim=0)
        mirna_weights = torch.norm(model.mirna_encoder[0].weight.data, dim=0)
        embedding_weights = torch.norm(model.embedding_encoder[0].weight.data, dim=0)
        
        all_weights = torch.cat([mutation_weights, mirna_weights, embedding_weights])
        max_weight = torch.max(all_weights)
        weights_scaled = (all_weights / max_weight * 100).cpu().numpy()
    
    n_features = len(feature_names)
    positions = np.arange(n_features) * 2
    
    for i, (feature, weight) in enumerate(zip(feature_names, weights_scaled)):
        plt.plot([0, weight], [positions[i], positions[i]], 'b-', linewidth=2)
        plt.text(-5, positions[i], feature, ha='right', va='center')
        
        ticks = np.linspace(0, weight, 5)
        plt.plot([ticks, ticks], 
                [positions[i]-0.1, positions[i]+0.1], 'k-')
    
    plt.ylim(min(positions)-2, max(positions)+2)
    plt.xlim(-10, 110)
    
    plt.plot([0, 100], [-1, -1], 'k-', linewidth=2)
    plt.text(50, -1.5, 'Points', ha='center')
    
    plt.title('Nomogram for Survival Prediction')
    plt.axis('off')
    
    return plt.gcf()

def plot_calibration_curves(model, data_loader, device, time_points=[365, 3*365, 5*365]):
    """Plot calibration curves for different time points"""
    plt.figure(figsize=(10, 8))
    
    model.eval()
    all_preds = []
    all_times = []
    all_events = []
    
    with torch.no_grad():
        for x_batch, time_batch, event_batch in data_loader:
            x_batch = [x.to(device) for x in x_batch]
            risk_pred = model(x_batch)
            all_preds.extend(risk_pred.squeeze().cpu().numpy())
            all_times.extend(time_batch.numpy())
            all_events.extend(event_batch.numpy())
    
    all_preds = np.array(all_preds)
    all_times = np.array(all_times)
    all_events = np.array(all_events)
    
    surv_probs = 1 / (1 + np.exp(all_preds))
    
    for time_point in time_points:
        actual = (all_times >= time_point) & (all_events == 1)
        if sum(actual) > 0:
            prob_true, prob_pred = calibration_curve(
                actual, surv_probs, n_bins=10, strategy='quantile'
            )
            plt.plot(prob_pred, prob_true, 
                    marker='o', 
                    label=f'{time_point/365:.1f} years')
    
    plt.plot([0, 1], [0, 1], 'k--', label='Ideal')
    plt.xlabel('Predicted Survival Probability')
    plt.ylabel('Actual Survival Probability')
    plt.title('Calibration Curves')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    return plt.gcf()

def calculate_metrics(model, data_loader, device):
    model.eval()
    all_preds = []
    all_times = []
    all_events = []
    
    with torch.no_grad():
        for x_batch, time_batch, event_batch in data_loader:
            x_batch = [x.to(device) for x in x_batch]
            risk_pred = model(x_batch)
            all_preds.extend(risk_pred.squeeze().cpu().numpy())
            all_times.extend(time_batch.numpy())
            all_events.extend(event_batch.numpy())
    
    all_preds = np.array(all_preds)
    all_times = np.array(all_times)
    all_events = np.array(all_events)
    
    c_index = concordance_index(all_times, -all_preds, all_events)
    
    times = np.array([365, 3*365, 5*365])
    structured_array = np.zeros(len(all_times), dtype=[('event', bool), ('time', float)])
    structured_array['event'] = all_events.astype(bool)
    structured_array['time'] = all_times
    
    try:
        _, mean_auc = cumulative_dynamic_auc(
            survival_train=structured_array,
            survival_test=structured_array,
            estimate=all_preds,
            times=times
        )
    except:
        mean_auc = 0.5
    
    return c_index, mean_auc

def train_fold(model, train_loader, val_loader, fold_num, n_epochs=100, patience=20):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.1)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='max', factor=0.5, patience=10, verbose=True
    )
    
    best_val_ci = 0
    best_val_auc = 0
    best_model = None
    patience_counter = 0
    
    for epoch in range(n_epochs):
        model.train()
        train_loss = 0
        
        for x_batch, time_batch, event_batch in train_loader:
            x_batch = [x.to(device) for x in x_batch]
            time_batch = time_batch.to(device)
            event_batch = event_batch.to(device)
            
            optimizer.zero_grad()
            risk_pred = model(x_batch)
            loss = cox_loss(risk_pred.squeeze(), time_batch, event_batch)
            
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            
            train_loss += loss.item()
        
        val_ci, val_auc = calculate_metrics(model, val_loader, device)
        scheduler.step(val_ci)
        
        if val_ci > best_val_ci:
            best_val_ci = val_ci
            best_val_auc = val_auc
            best_model = model.state_dict().copy()
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= patience:
            print(f'Fold {fold_num} - Early stopping at epoch {epoch+1}')
            print(f'Best Validation C-index: {best_val_ci:.4f}')
            print(f'Best Validation AUC: {best_val_auc:.4f}')
            break
        
        if (epoch + 1) % 10 == 0:
            print(f'Fold {fold_num} - Epoch {epoch+1}')
            print(f'Train Loss: {train_loss/len(train_loader):.4f}')
            print(f'Validation C-index: {val_ci:.4f}')
            print(f'Validation AUC: {val_auc:.4f}')
    
    return best_model, best_val_ci, best_val_auc

def evaluate_final_model(model, test_loader, device, output_dir=None):
    """Evaluate final model with nomogram and calibration curves"""
    model.eval()
    test_ci, test_auc = calculate_metrics(model, test_loader, device)
    
    print("\nFinal Test Set Results:")
    print(f"Test C-index: {test_ci:.4f}")
    print(f"Test AUC: {test_auc:.4f}")
    
    if output_dir:
        # Create feature names for nomogram
        feature_names = (
            [f'Mutation_{i}' for i in range(model.mutation_encoder[0].weight.shape[1])] +
            [f'miRNA_{i}' for i in range(model.mirna_encoder[0].weight.shape[1])] +
            [f'Embedding_{i}' for i in range(model.embedding_encoder[0].weight.shape[1])]
        )
        
        # Create and save nomogram
        nom_fig = create_nomogram(model, feature_names, device)
        nom_fig.savefig(f'{output_dir}/nomogram.png', bbox_inches='tight', dpi=300)
        plt.close(nom_fig)
        
        # Create and save calibration curves
        cal_fig = plot_calibration_curves(model, test_loader, device)
        cal_fig.savefig(f'{output_dir}/calibration_curves.png', bbox_inches='tight', dpi=300)
        plt.close(cal_fig)
    
    return test_ci, test_auc

def cross_validate_model(mutation_data, mirna_data, embedding_data, labels, n_splits=5):
    print(f"Starting {n_splits}-fold cross-validation")
    
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    cv_results = []
    best_overall_ci = 0
    best_overall_auc = 0
    best_fold = 0
    best_overall_model = None
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(mutation_data), 1):
        print(f"\nTraining Fold {fold}")
        
        train_dataset = MultiModalSurvivalDataset(
            mutation_data.iloc[train_idx], 
            mirna_data.iloc[train_idx], 
            embedding_data.iloc[train_idx], 
            labels.iloc[train_idx]
        )
        val_dataset = MultiModalSurvivalDataset(
            mutation_data.iloc[val_idx], 
            mirna_data.iloc[val_idx], 
            embedding_data.iloc[val_idx], 
            labels.iloc[val_idx]
        )
        
        train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
        
        input_dims = [
            mutation_data.shape[1],
            mirna_data.shape[1],
            embedding_data.shape[1]
        ]
        model = MultiModalDeepSurv(input_dims)
        
        best_model, best_ci, best_auc = train_fold(
            model, train_loader, val_loader, fold
        )
        
        cv_results.append({
            'fold': fold,
            'c_index': best_ci,
            'auc': best_auc
        })
        
        if best_ci > best_overall_ci:
            best_overall_ci = best_ci
            best_overall_auc = best_auc
            best_fold = fold
            best_overall_model = best_model
    
    ci_values = [result['c_index'] for result in cv_results]
    auc_values = [result['auc'] for result in cv_results]
    
    mean_ci = np.mean(ci_values)
    std_ci = np.std(ci_values)
    mean_auc = np.mean(auc_values)
    std_auc = np.std(auc_values)
    
    print(f"\nCross-validation summary:")
    print(f"Mean C-index: {mean_ci:.4f} ± {std_ci:.4f}")
    print(f"Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    print(f"\nBest Overall Results (Fold {best_fold}):")
    print(f"Best C-index: {best_overall_ci:.4f}")
    print(f"Best AUC: {best_overall_auc:.4f}")
    
    return {
        'mean_ci': mean_ci,
        'std_ci': std_ci,
        'mean_auc': mean_auc,
        'std_auc': std_auc,
        'best_model': best_overall_model,
        'best_ci': best_overall_ci,
        'best_auc': best_overall_auc,
        'cv_results': cv_results
    }

def main():
    # File paths
    mutation_file='/Users/stanleychen/git/Melanoma/data/pathway_maf_scores_filtered_above_cutoff.csv'
    mirna_file='/Users/stanleychen/git/Melanoma/apaper/miRNA.csv'
    embedding_file='/Users/stanleychen/git/Melanoma/final_data/final_embeddings.npy'
    id_mapping_file='/Users/stanleychen/git/Melanoma/final_data/final_identifiers.tsv'
    labels_file='/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv'
    output_dir='/Users/stanleychen/git/Melanoma/multimodaldeepsurv'
    
     # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Load and process data
        print("Loading and processing data...")
        mutations, mirna, embeddings, labels = process_and_align_data(
            mutation_file,
            mirna_file,
            embedding_file,
            id_mapping_file,
            labels_file
        )
        # Split into train+val and test sets
        train_val_idx, test_idx = train_test_split(
            np.arange(len(mutations)), 
            test_size=0.1, 
            random_state=42,
            stratify=labels['event']
        )
        
        # Create train+val and test sets
        mutations_train_val = mutations.iloc[train_val_idx]
        mirna_train_val = mirna.iloc[train_val_idx]
        embeddings_train_val = embeddings.iloc[train_val_idx]
        labels_train_val = labels.iloc[train_val_idx]
        
        mutations_test = mutations.iloc[test_idx]
        mirna_test = mirna.iloc[test_idx]
        embeddings_test = embeddings.iloc[test_idx]
        labels_test = labels.iloc[test_idx]
        
        print("\nDataset splits:")
        print(f"Train+Val set size: {len(mutations_train_val)}")
        print(f"Test set size: {len(mutations_test)}")
        
        # Perform cross-validation
        cv_results = cross_validate_model(
            mutations_train_val,
            mirna_train_val,
            embeddings_train_val,
            labels_train_val
        )
        
        # Create test dataset
        test_dataset = MultiModalSurvivalDataset(
            mutations_test,
            mirna_test,
            embeddings_test,
            labels_test
        )
        test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
        
        # Initialize and evaluate best model
        input_dims = [
            mutations.shape[1],
            mirna.shape[1],
            embeddings.shape[1]
        ]
        final_model = MultiModalDeepSurv(input_dims)
        final_model.load_state_dict(cv_results['best_model'])
        
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        final_model = final_model.to(device)
        
        # Evaluate with nomogram and calibration curves
        test_ci, test_auc = evaluate_final_model(
            final_model, 
            test_loader, 
            device,
            output_dir=output_dir
        )
        
        # Save results properly
        # Save CV results
        cv_results_df = pd.DataFrame([{
            'fold': result['fold'],
            'c_index': result['c_index'],
            'auc': result['auc']
        } for result in cv_results['cv_results']])
        cv_results_df.to_csv(f'{output_dir}/cv_results.csv', index=False)
        
        # Save summary results
        summary_results = pd.DataFrame({
            'metric': ['mean_ci', 'std_ci', 'mean_auc', 'std_auc', 'best_ci', 'best_auc', 'test_ci', 'test_auc'],
            'value': [
                cv_results['mean_ci'],
                cv_results['std_ci'],
                cv_results['mean_auc'],
                cv_results['std_auc'],
                cv_results['best_ci'],
                cv_results['best_auc'],
                test_ci,
                test_auc
            ]
        })
        summary_results.to_csv(f'{output_dir}/summary_results.csv', index=False)
        
        # Save data dimensions
        dimensions_df = pd.DataFrame({
            'modality': ['mutations', 'mirna', 'embeddings'],
            'samples': [mutations.shape[0], mirna.shape[0], embeddings.shape[0]],
            'features': [mutations.shape[1], mirna.shape[1], embeddings.shape[1]]
        })
        dimensions_df.to_csv(f'{output_dir}/data_dimensions.csv', index=False)
        
        print("\nResults saved to:", output_dir)
        return {
            'cv_results': cv_results,
            'test_results': {
                'c_index': test_ci,
                'auc': test_auc
            },
            'data_dimensions': {
                'mutations': mutations.shape,
                'mirna': mirna.shape,
                'embeddings': embeddings.shape
            }
        }
        
    except Exception as e:
        print(f"Error during execution: {str(e)}")
        raise

if __name__ == "__main__":
    main()