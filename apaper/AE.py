import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import warnings
warnings.filterwarnings('ignore')

class ExpressionDataset(Dataset):
    def __init__(self, data):
        self.data = torch.FloatTensor(data)
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        return self.data[idx]

class EarlyFusionAutoencoder(nn.Module):
    def __init__(self, input_dim):
        super(EarlyFusionAutoencoder, self).__init__()
        
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 500),
            nn.Tanh(),
            nn.Dropout(0.5),
            
            nn.Linear(500, 210),
            nn.Tanh(),
            nn.Dropout(0.5),
        )
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(210, 500),
            nn.Tanh(),
            nn.Dropout(0.5),
            
            nn.Linear(500, input_dim),
            nn.Tanh()
        )
        
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return encoded, decoded


def train_autoencoder(model, train_loader, device, epochs=50):
    # Initialize optimizer with L1 and L2 regularization
    optimizer = torch.optim.Adam(
        model.parameters(), 
        lr=0.001, 
        weight_decay=0.001  # L2 regularization
    )
    
    model.to(device)
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        
        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            
            encoded, decoded = model(batch)
            
            # MSE reconstruction loss
            mse_loss = F.mse_loss(decoded, batch)
            
            # L1 regularization
            l1_loss = 0
            for param in model.parameters():
                l1_loss += torch.sum(torch.abs(param))
            l1_loss *= 0.001  # L1 regularization strength
            
            # Total loss
            loss = mse_loss + l1_loss
            
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
        print(f'Epoch {epoch+1}, Loss: {total_loss/len(train_loader):.4f}')
    
    return model

def get_latent_space(model, data_loader, device):
    model.eval()
    all_encoded = []
    
    with torch.no_grad():
        for batch in data_loader:
            batch = batch.to(device)
            encoded, _ = model(batch)
            all_encoded.append(encoded.cpu().numpy())
    
    return np.vstack(all_encoded)

def calculate_survival_metrics(durations, events, cluster_labels, latent_features, patient_ids):
    """Calculate survival metrics for the clusters."""
    medians = []
    kmf = KaplanMeierFitter()

    for cluster in [0, 1]:
        mask = patient_ids.isin(latent_features[cluster_labels == cluster].index)
        kmf.fit(durations[mask], events[mask])
        medians.append(kmf.median_survival_time_)

        print(f"\nCluster {cluster + 1}:")
        print(f"Samples: {sum(mask)}")
        print(f"Median survival: {medians[-1]:.1f} days")
        print(f"Events (deaths): {sum(events[mask])}")
        print(f"Event rate: {sum(events[mask])/sum(mask):.2f}")

    high_risk_cluster = 0 if medians[0] < medians[1] else 1

    kmeans = KMeans(n_clusters=2, random_state=42)
    kmeans.fit(latent_features)
    distances = kmeans.transform(latent_features)
    risk_scores = distances[:, high_risk_cluster]

    c_index = concordance_index(durations, -risk_scores, events)

    high_risk_mask = cluster_labels == high_risk_cluster
    low_risk_mask = ~high_risk_mask
    p_value = logrank_test(
        durations[low_risk_mask], 
        durations[high_risk_mask],
        events[low_risk_mask],
        events[high_risk_mask]
    ).p_value

    print("\nSurvival Analysis:")
    print(f"High-risk cluster: Cluster {high_risk_cluster + 1}")
    print(f"C-index: {c_index:.3f}")
    print(f"Log-rank p-value: {p_value:.2e}")

    return c_index, p_value, high_risk_cluster

def process_and_align_data(expression_file, labels_file):
    """Load and align expression data with patient IDs and apply min-max scaling."""
    print("Loading and processing data...")

    expression_data = pd.read_csv(expression_file, index_col=0)
    print(f"Loaded expression data: {expression_data.shape}")
    
    # Print sample of raw expression data
    print("\nSample of raw expression data (first 5 genes, first 3 patients):")
    print(expression_data.iloc[:5, :3])

    labels = pd.read_csv(labels_file)
    labels = labels.drop_duplicates(subset='patient_barcode', keep='first')
    labels = labels.set_index('patient_barcode')
    
    # Print sample of labels
    print("\nSample of survival labels (first 5 patients):")
    print(labels[['survival_time_days', 'event', 'vital_status']].head())

    common_patients = sorted(list(set(expression_data.columns) & set(labels.index)))
    print(f"\nFound {len(common_patients)} common patients")

    if len(common_patients) == 0:
        raise ValueError("No common patients found between expression data and labels!")

    expression_data = expression_data[common_patients].T
    labels = labels.reindex(common_patients)

    scaler = MinMaxScaler()
    expression_scaled = scaler.fit_transform(expression_data)
    expression_data_scaled = pd.DataFrame(expression_scaled, index=expression_data.index, columns=expression_data.columns)
    
    # Print sample of processed data with corresponding labels
    print("\nSample of processed data with corresponding labels (first 3 patients, first 5 genes):")
    sample_df = pd.DataFrame({
        'Patient': expression_data_scaled.index[:3],
        'Survival_Days': labels['survival_time_days'][:3],
        'Event': labels['event'][:3],
        'Vital_Status': labels['vital_status'][:3]
    })
    print("\nPatient Information:")
    print(sample_df)
    
    print("\nCorresponding Expression Values (scaled):")
    print(expression_data_scaled.iloc[:3, :5])

    print(f"\nFinal expression data shape (scaled): {expression_data_scaled.shape}")
    
    # Print value ranges
    print("\nValue ranges:")
    print(f"Raw data - Min: {expression_data.values.min():.4f}, Max: {expression_data.values.max():.4f}")
    print(f"Scaled data - Min: {expression_data_scaled.values.min():.4f}, Max: {expression_data_scaled.values.max():.4f}")
    
    return expression_data_scaled, labels

def main():
    # File paths (modify as needed)
    expression_file = '/Users/stanleychen/git/Melanoma/apaper/RNA.csv'
    labels_file = '/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv'

    print("="*80)
    print("Starting data processing and analysis...")
    print("="*80)
    
    expression_data, labels = process_and_align_data(expression_file, labels_file)
    
    # Rest of the code remains the same...
    dataset = ExpressionDataset(expression_data.values)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"\nUsing device: {device}")
    
    model = EarlyFusionAutoencoder(expression_data.shape[1])
    model = train_autoencoder(model, dataloader, device, epochs=10)
    
    latent_features = get_latent_space(model, dataloader, device)
    latent_features = pd.DataFrame(latent_features, index=expression_data.index)

    kmeans = KMeans(n_clusters=2, random_state=42)
    cluster_labels = kmeans.fit_predict(latent_features.values)

    c_index, p_value, high_risk_cluster = calculate_survival_metrics(
        labels['survival_time_days'].values,
        labels['event'].values,
        cluster_labels,
        latent_features,
        expression_data.index
    )
    
    return {
        'c_index': c_index,
        'p_value': p_value,
        'cluster_labels': cluster_labels,
        'latent_features': latent_features,
        'model': model
    }

if __name__ == "__main__":
    results = main()