import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from sksurv.metrics import cumulative_dynamic_auc
import warnings
warnings.filterwarnings('ignore')

class DualModalDeepSurv(nn.Module):
    """
    Modified DualModal architecture to more closely match successful single modality version
    """
    def __init__(self, input_dims, drop=0.3):
        super(DualModalDeepSurv, self).__init__()
        
        # Individual modality networks following original DeepSurv architecture
        self.mirna_net = nn.Sequential(
            nn.Linear(input_dims[0], 64),
            nn.ReLU(),
            nn.Dropout(drop),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(drop)
        )
        
        self.embedding_net = nn.Sequential(
            nn.Linear(input_dims[1], 64),
            nn.ReLU(),
            nn.Dropout(drop),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(drop)
        )
        
        # Final fusion layer
        self.fusion_net = nn.Linear(64, 1)  # 32 + 32 = 64 combined features
        
        # Initialize weights according to paper
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                nn.init.zeros_(m.bias)
    
    def forward(self, x):
        mirna_data, embedding_data = x
        # Process each modality
        mirna_features = self.mirna_net(mirna_data)
        embedding_features = self.embedding_net(embedding_data)
        # Concatenate and final prediction
        combined = torch.cat([mirna_features, embedding_features], dim=1)
        return self.fusion_net(combined)

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

def cox_loss(risk_pred, y_time, y_event):
    """
    Calculate the Cox proportional hazards loss
    """
    # Sort by survival time in descending order
    order = torch.argsort(y_time, descending=True)
    risk_pred = risk_pred[order]
    y_event = y_event[order]
    
    # Calculate log risk
    log_risk = torch.log(torch.cumsum(torch.exp(risk_pred), dim=0))
    
    # Calculate partial likelihood
    likelihood = risk_pred - log_risk
    
    # Apply censoring mask
    partial_likelihood = likelihood * y_event
    
    # Return negative average partial likelihood
    return -torch.mean(partial_likelihood)

def calculate_time_dependent_auc(times, events, scores):
    """
    Calculate time-dependent AUC for survival data
    """
    # Convert to structured array
    structured_array = np.zeros(len(times), dtype=[('event', bool), ('time', float)])
    structured_array['event'] = events.astype(bool)
    structured_array['time'] = times
    
    # Define meaningful time points
    time_points = np.array([365, 3*365, 5*365])  # 1, 3, and 5 years
    
    # For each time point, ensure we have enough events
    valid_times = []
    
    print("\nTime-dependent AUC analysis:")
    for t in time_points:
        # Count events up to this time
        events_before_t = np.sum((times <= t) & events.astype(bool))
        at_risk_at_t = np.sum(times > t)
        
        print(f"\nAt {t/365:.1f} years:")
        print(f"- Events before: {events_before_t}")
        print(f"- At risk: {at_risk_at_t}")
        
        if events_before_t >= 5 and at_risk_at_t >= 5:  # Minimum threshold
            valid_times.append(t)
    
    if len(valid_times) == 0:
        print("Warning: No valid time points found with sufficient events")
        return 0.5, np.array([0.5])
    
    try:
        with np.errstate(divide='ignore', invalid='ignore'):
            auc, mean_auc = cumulative_dynamic_auc(
                survival_train=structured_array,
                survival_test=structured_array,
                estimate=scores,
                times=valid_times
            )
            
            # Handle NaN values
            auc = np.nan_to_num(auc, nan=0.5)
            mean_auc = np.nan_to_num(mean_auc, nan=0.5)
            
            # Print AUC for each time point
            for t, auc_t in zip(valid_times, auc):
                print(f"\nAUC at {t/365:.1f} years: {auc_t:.3f}")
            
    except Exception as e:
        print(f"Warning: Could not calculate time-dependent AUC: {str(e)}")
        return 0.5, np.array([0.5])
    
    return mean_auc, auc

def process_and_align_data(mirna_file, embedding_file, id_mapping_file, labels_file):
    """
    Align miRNA and RNA-seq data using patient IDs
    """
    print("Loading and processing data...")
    
    # Load ID mappings
    id_mappings = pd.read_csv(id_mapping_file, sep='\t', names=['uuid', 'patient_id'])
    id_mappings = id_mappings.drop_duplicates(subset='patient_id', keep='first')
    uuid_to_patient = dict(zip(id_mappings['uuid'], id_mappings['patient_id']))
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
    
    # Load embeddings from .npy file
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
    
    # Find common patients across modalities
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

def train_model(model, train_loader, val_loader, n_epochs=500, lr=0.01, patience=50):
    """
    Train the DualModal DeepSurv model with matched hyperparameters
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    optimizer = torch.optim.SGD(model.parameters(), lr=lr, weight_decay=0.1)
    
    best_ci = 0
    best_auc = 0
    best_model = None
    patience_counter = 0
    
    for epoch in range(n_epochs):
        # Training
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
            optimizer.step()
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_preds = []
        val_times = []
        val_events = []
        
        with torch.no_grad():
            for x_batch, time_batch, event_batch in val_loader:
                x_batch = [x.to(device) for x in x_batch]
                risk_pred = model(x_batch)
                val_preds.extend(risk_pred.squeeze().cpu().numpy())
                val_times.extend(time_batch.numpy())
                val_events.extend(event_batch.numpy())
        
        val_preds = np.array(val_preds)
        val_times = np.array(val_times)
        val_events = np.array(val_events)
        
        # Calculate metrics
        try:
            val_ci = concordance_index(val_times, -val_preds, val_events)
        except:
            val_ci = 0.5
            
        try:
            mean_auc, auc_times = calculate_time_dependent_auc(
                val_times,
                val_events,
                val_preds
            )
        except:
            mean_auc = 0.5
            auc_times = np.array([0.5])
        
        # Update best model
        if val_ci > best_ci:
            best_ci = val_ci
            best_auc = mean_auc
            best_model = model.state_dict().copy()
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= patience:
            print(f'Early stopping at epoch {epoch+1}')
            break
        
        if (epoch + 1) % 50 == 0:
            print(f'Epoch {epoch+1}/{n_epochs}')
            print(f'Train Loss: {train_loss/len(train_loader):.4f}')
            print(f'Validation C-index: {val_ci:.4f}')
            print(f'Validation Mean AUC: {mean_auc:.4f}')
    
    return best_model, best_ci, best_auc

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
    
    # Create dataset
    dataset = DualModalSurvivalDataset(mirna, embeddings, labels)
    
    # Split data with same random seed as single modality
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(
        dataset, 
        [train_size, val_size],
        generator=torch.Generator().manual_seed(42)
    )
    
    # Create dataloaders with matched batch size
    train_loader = DataLoader(train_dataset, batch_size=100, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=100, shuffle=False)
    
    # Initialize model
    input_dims = [mirna.shape[1], embeddings.shape[1]]
    model = DualModalDeepSurv(input_dims)
    
    # Train model
    print("\nStarting training...")
    best_model, best_ci, best_auc = train_model(model, train_loader, val_loader)
    
    print("\nTraining completed!")
    print(f"Best C-index: {best_ci:.4f}")
    print(f"Best AUC: {best_auc:.4f}")
    
    # Save model and results
    results = {
        'model_state_dict': best_model,
        'c_index': best_ci,
        'auc': best_auc,
        'input_dims': input_dims,
        'feature_dims': {
            'mirna': mirna.shape[1],
            'embeddings': embeddings.shape[1]
        }
    }
    
    torch.save(results, 'dual_modal_deepsurv_model.pt')
    print("\nModel saved as 'dual_modal_deepsurv_model.pt'")
    
    return results

if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"Error in main execution: {str(e)}")
        raise