import torch
import torch.nn as nn
import numpy as np
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from lifelines.utils import concordance_index
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sksurv.metrics import cumulative_dynamic_auc
from sksurv.util import Surv

class DeepSurv(nn.Module):
    """
    Original DeepSurv architecture from Katzman et al. 2016
    """
    def __init__(self, input_dim, drop=0.3):
        super(DeepSurv, self).__init__()
        # Architecture from the original paper
        self.net = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Dropout(drop),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(drop),
            nn.Linear(32, 1)
        )
        
        # Initialize weights according to paper
        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.kaiming_normal_(m.weight)
                nn.init.zeros_(m.bias)
    
    def forward(self, x):
        return self.net(x)

def cox_loss(risk_pred, y_time, y_event):
    """
    True Cox partial likelihood loss from the DeepSurv paper
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
    Calculate time-dependent AUC for survival data with specific handling for censored data
    
    Parameters:
    -----------
    times : array-like
        Survival times in days
    events : array-like
        Event indicators (1 for death, 0 for censored/alive)
    scores : array-like
        Predicted risk scores from the model
    """
    # Convert to structured array
    structured_array = np.zeros(len(times), dtype=[('event', bool), ('time', float)])
    structured_array['event'] = events.astype(bool)
    structured_array['time'] = times
    
    # Define meaningful time points based on your data distribution
    # Using 1-year, 3-year, and 5-year survival as key points
    time_points = np.array([365, 3*365, 5*365])  # 1, 3, and 5 years
    
    # For each time point, ensure we have enough events
    valid_times = []
    aucs = []
    
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

class SurvivalDataset(Dataset):
    def __init__(self, X, time, event):
        """
        Initialize dataset with proper conversion from pandas/numpy to torch
        """
        # Convert to numpy first if needed
        if isinstance(X, pd.DataFrame) or isinstance(X, pd.Series):
            X = X.values
        if isinstance(time, pd.Series):
            time = time.values
        if isinstance(event, pd.Series):
            event = event.values
            
        # Convert to torch tensors
        self.X = torch.FloatTensor(X)
        self.time = torch.FloatTensor(time)
        self.event = torch.FloatTensor(event)
        
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.time[idx], self.event[idx]

def train_deepsurv(model, train_loader, val_loader, n_epochs=500, lr=0.01, patience=50):
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
        for X_batch, time_batch, event_batch in train_loader:
            X_batch = X_batch.to(device)
            time_batch = time_batch.to(device)
            event_batch = event_batch.to(device)
            
            optimizer.zero_grad()
            risk_pred = model(X_batch)
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
            for X_batch, time_batch, event_batch in val_loader:
                X_batch = X_batch.to(device)
                risk_pred = model(X_batch)
                val_preds.extend(risk_pred.squeeze().cpu().numpy())
                val_times.extend(time_batch.numpy())
                val_events.extend(event_batch.numpy())
        
        val_preds = np.array(val_preds)
        val_times = np.array(val_times)
        val_events = np.array(val_events)
        
        # Calculate metrics
        try:
            val_ci = concordance_index(
                val_times,
                -val_preds,
                val_events
            )
        except:
            val_ci = 0.5  # Default value if calculation fails
            
        try:
            mean_auc, auc_times = calculate_time_dependent_auc(
                val_times,
                val_events,
                val_preds
            )
        except:
            mean_auc = 0.5  # Default value if calculation fails
            auc_times = np.array([0.5])
        
        # Update best model based on C-index only (more reliable for survival)
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
    # Load feature data
    features_df = pd.read_csv('/Users/stanleychen/git/Melanoma/maf_pipeline/patient_maf_matrix_500.csv', index_col=0)

    # Transpose the data so that patient identifiers are rows
    # features_df = features_df.T

    # Reset the index to make patient identifiers a column
    features_df = features_df.reset_index()

    # Rename the identifier column to match the label file
    features_df.rename(columns={'index': 'patient_barcode'}, inplace=True)

    # Load survival data
    survival_df = pd.read_csv('/Users/stanleychen/git/Melanoma/data/patient_survival_analysis_filtered_keep_negs.csv')

    # Merge data on 'PatientID'
    data = pd.merge(features_df, survival_df, on='patient_barcode')
    data['vital_status'] = data['vital_status'].map({'Dead': 1, 'Alive': 0})

    # Extract features, survival time, and event status
    #X = data.drop(columns=['patient_barcode', 'survival_time_days', 'vital_status', 'survival_time_source','tumor_stage','grade'])
    X = data.drop(columns=['patient_barcode', 'survival_time_days', 'vital_status', 'survival_time_source','tumor_stage'])
    time = data['survival_time_days']
    event = data['vital_status']
    
    print(f"Dataset shape - Features: {X.shape}, Time: {len(time)}, Events: {len(event)}")
    
    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Split data
    print("Splitting data into train/val sets...")
    X_train, X_val, time_train, time_val, event_train, event_val = train_test_split(
        X_scaled, time, event, test_size=0.2, random_state=42
    )
    
    # Create datasets
    train_dataset = SurvivalDataset(X_train, time_train, event_train)
    val_dataset = SurvivalDataset(X_val, time_val, event_val)
    
    # Create dataloaders
    train_loader = DataLoader(train_dataset, batch_size=100, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=100, shuffle=False)
    
    # Initialize model
    print("Initializing model...")
    model = DeepSurv(input_dim=X.shape[1])
    
    # Train model
    print("Starting training...")
    best_model, best_ci, best_auc = train_deepsurv(model, train_loader, val_loader)
    print(f'\nFinal Results:')
    print(f'Best validation C-index: {best_ci:.4f}')
    print(f'Best validation AUC: {best_auc:.4f}')
    
    # Save model and metrics
    torch.save({
        'model_state_dict': best_model,
        'c_index': best_ci,
        'auc': best_auc
    }, 'deepsurv_model.pt')
    print("Model saved as 'deepsurv_model.pt'")

if __name__ == '__main__':
    main()