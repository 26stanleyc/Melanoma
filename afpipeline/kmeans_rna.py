import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import warnings
warnings.filterwarnings('ignore')

class RNAseqDataset(Dataset):
    def __init__(self, data):
        self.data = torch.FloatTensor(data)
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        return self.data[idx]

class Autoencoder(nn.Module):
    def __init__(self, input_dim):
        super(Autoencoder, self).__init__()
        
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(0.4),
            
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.4),
            
            nn.Linear(256, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            
            nn.Linear(64, 2)  # 2D latent space
        )
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(2, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            
            nn.Linear(64, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            
            nn.Linear(256, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            
            nn.Linear(512, input_dim)
        )
        
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return encoded, decoded

def process_and_align_data(embedding_file, id_mapping_file, labels_file):
    """Load and align RNA-seq data with patient IDs"""
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
    
    # Filter embeddings to matched patients
    embeddings = embeddings[matched_indices]
    embeddings_df = pd.DataFrame(
        embeddings,
        index=matched_patients,
        columns=[f'emb_{i}' for i in range(embeddings.shape[1])]
    )
    print(f"Processed embeddings shape: {embeddings_df.shape}")
    
    # Keep only patients with labels
    embeddings_df = embeddings_df[embeddings_df.index.isin(valid_patients)]
    labels = labels[labels['patient_barcode'].isin(embeddings_df.index)]
    
    # Convert survival times and events
    labels['event'] = (labels['vital_status'] == 'Dead').astype(int)
    
    # Sort dataframes by patient ID
    common_patients = sorted(list(set(embeddings_df.index) & set(labels['patient_barcode'])))
    embeddings_df = embeddings_df.reindex(index=common_patients)
    labels = labels.set_index('patient_barcode').reindex(index=common_patients)
    
    print(f"Final number of samples: {len(common_patients)}")
    
    return embeddings_df, labels

def train_autoencoder(model, train_loader, device, epochs=100):
    """Train autoencoder model"""
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    model.to(device)
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        
        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            
            encoded, decoded = model(batch)
            loss = F.mse_loss(decoded, batch)
            
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
        if (epoch + 1) % 10 == 0:
            print(f'Epoch {epoch+1}, Loss: {total_loss/len(train_loader):.4f}')
    
    return model

def get_latent_space(model, data_loader, device):
    """Get latent space representations"""
    model.eval()
    all_encoded = []
    
    with torch.no_grad():
        for batch in data_loader:
            batch = batch.to(device)
            encoded, _ = model(batch)
            all_encoded.append(encoded.cpu().numpy())
    
    return np.vstack(all_encoded)

def calculate_survival_metrics(durations, events, cluster_labels, latent_features):
    """Calculate survival metrics with proper risk assignment"""
    # Calculate median survival for each cluster
    medians = []
    kmf = KaplanMeierFitter()
    for cluster in [0, 1]:
        mask = cluster_labels == cluster
        kmf.fit(durations[mask], events[mask])
        medians.append(kmf.median_survival_time_)
        
        print(f"\nCluster {cluster + 1}:")
        print(f"Samples: {sum(mask)}")
        print(f"Median survival: {medians[-1]:.1f} days")
        print(f"Events (deaths): {sum(events[mask])}")
        print(f"Event rate: {sum(events[mask])/sum(mask):.2f}")
    
    # Identify high risk cluster (lower survival = higher risk)
    high_risk_cluster = 0 if medians[0] < medians[1] else 1
    
    # Get distances to cluster centers for risk scores
    kmeans = KMeans(n_clusters=2, random_state=42)
    kmeans.fit(latent_features)
    distances = kmeans.transform(latent_features)
    risk_scores = distances[:, high_risk_cluster]
    
    # Calculate C-index
    c_index = concordance_index(durations, risk_scores, events)
    
    # Calculate log-rank p-value
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

def plot_latent_space(latent_features, cluster_labels, high_risk_cluster, output_path):
    """Plot latent space with risk groups"""
    plt.figure(figsize=(10, 8))
    
    colors = ['#FF9999', '#66CCCC']  # High-risk: red, Low-risk: blue
    
    for i in range(2):
        mask = cluster_labels == i
        points = latent_features[mask]
        risk_label = 'High-risk' if i == high_risk_cluster else 'Low-risk'
        color = colors[0] if i == high_risk_cluster else colors[1]
        
        plt.scatter(points[:, 0], points[:, 1],
                   c=[color], alpha=0.6, s=50,
                   marker='^' if i == 1 else 'o',
                   label=f'{risk_label} (n={sum(mask)})')
        
        # Create convex hull
        if len(points) >= 3:
            hull = ConvexHull(points)
            hull_vertices = points[hull.vertices]
            hull_vertices = np.vstack([hull_vertices, hull_vertices[0]])
            plt.fill(hull_vertices[:, 0], hull_vertices[:, 1],
                    alpha=0.1, c=color)
    
    plt.grid(True, alpha=0.3)
    plt.legend(title='Risk Groups')
    plt.title('RNA-seq Based Clustering')
    plt.xlabel('Latent Dimension 1')
    plt.ylabel('Latent Dimension 2')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_km_curves(durations, events, cluster_labels, high_risk_cluster, output_path):
    """Create Kaplan-Meier survival plots for risk groups"""
    plt.figure(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    colors = {'Low-risk': '#66CCCC', 'High-risk': '#FF9999'}
    
    for i in range(2):
        mask = cluster_labels == i
        risk_label = 'High-risk' if i == high_risk_cluster else 'Low-risk'
        
        kmf.fit(
            durations[mask],
            events[mask],
            label=f"{risk_label} (n={sum(mask)})"
        )
        kmf.plot(
            ci_show=True,
            color=colors[risk_label],
            alpha=0.7
        )
    
    plt.title('Survival Analysis of SKCM Risk Groups')
    plt.xlabel('Time (days)')
    plt.ylabel('Survival probability')
    plt.grid(True, alpha=0.3)
    
    # Add log-rank p-value
    high_risk_mask = cluster_labels == high_risk_cluster
    low_risk_mask = ~high_risk_mask
    p_value = logrank_test(
        durations[low_risk_mask], 
        durations[high_risk_mask],
        events[low_risk_mask],
        events[high_risk_mask]
    ).p_value
    
    plt.text(0.05, 0.05, f'Log-rank P = {p_value:.2e}', 
             transform=plt.gca().transAxes)
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()

def main():
    # File paths
    mutation_file = '/Users/stanleychen/git/Melanoma/data/pathway_maf_scores_filtered_above_cutoff.csv'
    mirna_file='/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA_after_IGCBA.csv'
    embedding_file='/Users/stanleychen/git/Melanoma/final_data/final_embeddings.npy'
    id_mapping_file='/Users/stanleychen/git/Melanoma/final_data/final_identifiers.tsv'
    labels_file='/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv'
    
    # Load and process data
    embedding_data, labels = process_and_align_data(
        embedding_file, id_mapping_file, labels_file
    )
    
    # Scale the data
    scaler = StandardScaler()
    embedding_scaled = scaler.fit_transform(embedding_data)
    
    # Create dataset and dataloader
    dataset = RNAseqDataset(embedding_scaled)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    # Device configuration
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"\nUsing device: {device}")
    
    # Train Autoencoder
    print("\nTraining Autoencoder...")
    model = Autoencoder(embedding_scaled.shape[1])
    model = train_autoencoder(model, dataloader, device)
    
    # Get latent representations
    latent_features = get_latent_space(model, dataloader, device)
    
    # Perform clustering
    print("\nPerforming clustering...")
    kmeans = KMeans(n_clusters=2, random_state=42)
    cluster_labels = kmeans.fit_predict(latent_features)
    
    # Calculate survival metrics
    c_index, p_value, high_risk_cluster = calculate_survival_metrics(
        labels['survival_time_days'].values,
        labels['event'].values,
        cluster_labels,
        latent_features
    )
    
    # Plot results
    plot_latent_space(
        latent_features,
        cluster_labels,
        high_risk_cluster,
        'rnaseq_clusters.png'
    )
    
    plot_km_curves(
        labels['survival_time_days'].values,
        labels['event'].values,
        cluster_labels,
        high_risk_cluster,
        'rnaseq_survival.png'
    )
    
    # Save cluster assignments
    cluster_df = pd.DataFrame({
        'patient_id': labels.index,
        'cluster': cluster_labels,
        'risk_group': ['High-risk' if c == high_risk_cluster else 'Low-risk' 
                      for c in cluster_labels],
        'survival_time': labels['survival_time_days'],
        'event': labels['event']
    })
    cluster_df.to_csv('rnaseq_clusters.csv')
    
    return {
        'c_index': c_index,
        'p_value': p_value,
        'cluster_labels': cluster_labels,
        'latent_features': latent_features
    }

if __name__ == "__main__":
    results = main()