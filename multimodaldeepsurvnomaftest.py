import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, train_test_split
import warnings
warnings.filterwarnings('ignore')

class DualModalSurvivalDataset(Dataset):
    def __init__(self, mirna_data, embedding_data, labels):
        """
        Initialize dataset with miRNA and RNA-seq embedding data
        """
        self.mirna_data = torch.FloatTensor(mirna_data.values)
        self.embedding_data = torch.FloatTensor(embedding_data.values)
        self.times = torch.FloatTensor(labels['survival_time_days'].values)
        self.events = torch.FloatTensor(labels['event'].values)
        
    def __len__(self):
        return len(self.times)
    
    def __getitem__(self, idx):
        return (self.mirna_data[idx], 
                self.embedding_data[idx]), \
               self.times[idx], \
               self.events[idx]

class SimpleDualModalDeepSurv(nn.Module):
    def __init__(self, input_dims, hidden_dim=32):
        super(SimpleDualModalDeepSurv, self).__init__()
        
        self.mirna_encoder = nn.Sequential(
            nn.Linear(input_dims[0], hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.4)
        )
        
        self.embedding_encoder = nn.Sequential(
            nn.Linear(input_dims[1], hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.4)
        )
        
        self.fusion_net = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(hidden_dim, 1)
        )
        
        self.apply(self._init_weights)
    
    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.normal_(m.weight, std=0.01)
            nn.init.zeros_(m.bias)
    
    def forward(self, x):
        mirna_data, embedding_data = x
        mirna_features = self.mirna_encoder(mirna_data)
        embedding_features = self.embedding_encoder(embedding_data)
        fused = torch.cat([mirna_features, embedding_features], dim=1)
        return self.fusion_net(fused)

def cox_loss(risk_pred, y_time, y_event):
    """
    Calculate the Cox proportional hazards loss
    """
    order = torch.argsort(y_time, descending=True)
    risk_pred = risk_pred[order]
    y_event = y_event[order]
    
    log_risk = torch.log(torch.cumsum(torch.exp(risk_pred), dim=0))
    likelihood = risk_pred - log_risk
    uncensored_likelihood = likelihood * y_event
    
    return -torch.mean(uncensored_likelihood)

def evaluate_model(model, data_loader, device):
    """
    Evaluate model on given dataset
    """
    model.eval()
    all_preds, all_times, all_events = [], [], []
    
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
    return c_index

def train_fold(model, train_loader, val_loader, fold_num):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, weight_decay=0.1)
    
    best_val_ci = 0
    best_train_ci = 0
    best_model = None
    patience = 20
    patience_counter = 0
    
    for epoch in range(100):
        # Training
        model.train()
        train_loss = 0
        train_preds, train_times, train_events = [], [], []
        
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
            train_preds.extend(risk_pred.squeeze().detach().cpu().numpy())
            train_times.extend(time_batch.cpu().numpy())
            train_events.extend(event_batch.cpu().numpy())
        
        # Calculate metrics
        train_ci = concordance_index(
            np.array(train_times), 
            -np.array(train_preds), 
            np.array(train_events)
        )
        
        # Validation
        val_ci = evaluate_model(model, val_loader, device)
        
        if val_ci > best_val_ci:
            best_val_ci = val_ci
            best_train_ci = train_ci
            best_model = model.state_dict().copy()
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= patience:
            print(f'Fold {fold_num} - Early stopping at epoch {epoch+1}')
            print(f'Best Training C-index: {best_train_ci:.4f}')
            print(f'Best Validation C-index: {best_val_ci:.4f}')
            break
        
        if (epoch + 1) % 10 == 0:
            print(f'Fold {fold_num} - Epoch {epoch+1}')
            print(f'Train Loss: {train_loss/len(train_loader):.4f}')
            print(f'Training C-index: {train_ci:.4f}')
            print(f'Validation C-index: {val_ci:.4f}')
    
    return best_model, best_val_ci, best_train_ci

def process_and_align_data(mirna_file, embedding_file, id_mapping_file, labels_file):
    """
    Align miRNA and RNA-seq data using patient IDs
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
    
    # Load miRNA data
    mirna = pd.read_csv(mirna_file, index_col=0)
    mirna = mirna[~mirna.index.duplicated(keep='first')]
    mirna = mirna[mirna.index.isin(valid_patients)]
    print(f"Loaded miRNA data: {mirna.shape}")
    
    # Load embeddings
    embeddings = np.load(embedding_file)
    print(f"Loaded embeddings with shape: {embeddings.shape}")
    
    # Match embeddings with patient IDs
    matched_indices = []
    matched_patients = []
    patient_to_idx = {}
    
    print("\nMatching embeddings to patient IDs...")
    for idx, (_, row) in enumerate(id_mappings.iterrows()):
        patient_id = row['patient_id']
        if patient_id in valid_patients and patient_id not in patient_to_idx:
            matched_indices.append(idx)
            matched_patients.append(patient_id)
            patient_to_idx[patient_id] = len(matched_indices) - 1
    
    if not matched_indices:
        raise ValueError("No matching patients found between embeddings and labels!")
    
    # Filter embeddings to matched patients
    embeddings = embeddings[matched_indices]
    embeddings_df = pd.DataFrame(
        embeddings,
        index=matched_patients,
        columns=[f'emb_{i}' for i in range(embeddings.shape[1])]
    )
    print(f"Processed embeddings shape: {embeddings_df.shape}")
    
    # Find common patients
    common_patients = set(mirna.index) & set(embeddings_df.index)
    print(f"Found {len(common_patients)} common patients across modalities")
    
    if len(common_patients) == 0:
        raise ValueError("No common patients found across modalities!")
    
    # Filter all data to common patients
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
    
    # Sort all dataframes by patient ID for consistency
    common_patients = sorted(list(common_patients))
    mirna_scaled = mirna_scaled.reindex(index=common_patients)
    embeddings_scaled = embeddings_scaled.reindex(index=common_patients)
    labels = labels.set_index('patient_barcode').reindex(index=common_patients)
    
    print("\nFinal dataset sizes after alignment:")
    print(f"miRNA: {mirna_scaled.shape}")
    print(f"Embeddings: {embeddings_scaled.shape}")
    print(f"Labels: {labels.shape}")
    
    return mirna_scaled, embeddings_scaled, labels

def cross_validate_model(mirna, embeddings, labels, n_splits=5):
    print(f"Starting {n_splits}-fold cross-validation")
    
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    cv_results = []
    train_results = []
    best_overall_ci = 0
    best_overall_train_ci = 0
    best_fold = 0
    best_overall_model = None
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(mirna.index), 1):
        print(f"\nTraining Fold {fold}")
        
        # Create fold-specific datasets
        train_dataset = DualModalSurvivalDataset(
            mirna.iloc[train_idx], 
            embeddings.iloc[train_idx], 
            labels.iloc[train_idx]
        )
        val_dataset = DualModalSurvivalDataset(
            mirna.iloc[val_idx], 
            embeddings.iloc[val_idx], 
            labels.iloc[val_idx]
        )
        
        # Create dataloaders
        train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
        
        # Initialize and train model
        input_dims = [mirna.shape[1], embeddings.shape[1]]
        model = SimpleDualModalDeepSurv(input_dims)
        
        best_model, best_ci, best_train_ci = train_fold(model, train_loader, val_loader, fold)
        cv_results.append(best_ci)
        train_results.append(best_train_ci)
        
        if best_ci > best_overall_ci:
            best_overall_ci = best_ci
            best_overall_train_ci = best_train_ci
            best_fold = fold
            best_overall_model = best_model
        
        print(f"Fold {fold} Best C-index: {best_ci:.4f}")
    
    mean_ci = np.mean(cv_results)
    std_ci = np.std(cv_results)
    mean_train_ci = np.mean(train_results)
    std_train_ci = np.std(train_results)
    
    print(f"\nCross-validation results:")
    print(f"Mean Training C-index: {mean_train_ci:.4f} ± {std_train_ci:.4f}")
    print(f"Mean Validation C-index: {mean_ci:.4f} ± {std_ci:.4f}")
    print(f"\nBest Overall Results (Fold {best_fold}):")
    print(f"Best Training C-index: {best_overall_train_ci:.4f}")
    print(f"Best Validation C-index: {best_overall_ci:.4f}")
    
    return mean_ci, std_ci, cv_results, best_overall_ci, best_overall_train_ci, best_overall_model

def main():
    # File paths
    mirna_file='/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA_after_IGCBA.csv'
    embedding_file='/Users/stanleychen/git/Melanoma/final_data/final_embeddings.npy'
    id_mapping_file='/Users/stanleychen/git/Melanoma/final_data/final_identifiers.tsv'
    labels_file='/Users/stanleychen/git/Melanoma/data/patient_survival_analysis_filtered_keep_negs.csv'
    
    # Load and align data
    try:
        mirna, embeddings, labels = process_and_align_data(
            mirna_file,
            embedding_file,
            id_mapping_file,
            labels_file
        )
    except Exception as e:
        print(f"Error during data processing: {str(e)}")
        raise
    
    # Split into train+val and test sets (90/10 split)
    train_val_idx, test_idx = train_test_split(
        np.arange(len(mirna)), 
        test_size=0.1, 
        random_state=42,
        stratify=labels['event']  # Stratify by event status
    )
    
    # Create DataFrames for train+val and test sets
    mirna_train_val = mirna.iloc[train_val_idx]
    embeddings_train_val = embeddings.iloc[train_val_idx]
    labels_train_val = labels.iloc[train_val_idx]
    
    mirna_test = mirna.iloc[test_idx]
    embeddings_test = embeddings.iloc[test_idx]
    labels_test = labels.iloc[test_idx]
    
    print("\nDataset splits:")
    print(f"Train+Val set size: {len(mirna_train_val)}")
    print(f"Test set size: {len(mirna_test)}")
    
    # Perform cross-validation on train+val set
    mean_ci, std_ci, cv_results, best_cv_ci, best_train_ci, best_model = cross_validate_model(
        mirna_train_val, 
        embeddings_train_val, 
        labels_train_val
    )
    
    # Create test dataset and loader
    test_dataset = DualModalSurvivalDataset(mirna_test, embeddings_test, labels_test)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    
    # Initialize best model for final evaluation
    input_dims = [mirna.shape[1], embeddings.shape[1]]
    final_model = SimpleDualModalDeepSurv(input_dims)
    final_model.load_state_dict(best_model)
    
    # Evaluate on test set
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    final_model = final_model.to(device)
    test_ci = evaluate_model(final_model, test_loader, device)
    
    print("\nFinal Results Summary:")
    print(f"Cross-Validation Results:")
    print(f"  Individual fold results: {[f'{ci:.4f}' for ci in cv_results]}")
    print(f"  Mean CV C-index: {mean_ci:.4f} ± {std_ci:.4f}")
    print(f"Best Model Performance:")
    print(f"  Best Training C-index: {best_train_ci:.4f}")
    print(f"  Best CV C-index: {best_cv_ci:.4f}")
    print(f"Test Set Performance:")
    print(f"  Test C-index: {test_ci:.4f}")

if __name__ == "__main__":
    main()