import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import ConvexHull
import warnings
warnings.filterwarnings('ignore')

class MultiOmicsDataset(Dataset):
    def __init__(self, mutation_data, mirna_data, embedding_data):
        self.mutation_data = torch.FloatTensor(mutation_data)
        self.mirna_data = torch.FloatTensor(mirna_data)
        self.embedding_data = torch.FloatTensor(embedding_data)
        
    def __len__(self):
        return len(self.mutation_data)
    
    def __getitem__(self, idx):
        return (self.mutation_data[idx], 
                self.mirna_data[idx], 
                self.embedding_data[idx])

class EarlyFusionAE(nn.Module):
    def __init__(self, input_dims):
        super(EarlyFusionAE, self).__init__()
        
        self.total_dim = sum(input_dims)
        
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(self.total_dim, 512),
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
            
            nn.Linear(512, self.total_dim)
        )
        
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return encoded, decoded

class LateFusionAE(nn.Module):
    def __init__(self, input_dims):
        super(LateFusionAE, self).__init__()
        
        mutation_dim, mirna_dim, embedding_dim = input_dims
        
        # Separate encoders for each modality
        self.mutation_encoder = nn.Sequential(
            nn.Linear(mutation_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Linear(256, 64)
        )
        
        self.mirna_encoder = nn.Sequential(
            nn.Linear(mirna_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Linear(256, 64)
        )
        
        self.embedding_encoder = nn.Sequential(
            nn.Linear(embedding_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Linear(256, 64)
        )
        
        # Fusion layer
        self.fusion = nn.Sequential(
            nn.Linear(64 * 3, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Linear(64, 2)  # 2D latent space
        )
        
        # Decoders
        self.mutation_decoder = nn.Sequential(
            nn.Linear(2, 64),
            nn.ReLU(),
            nn.Linear(64, 256),
            nn.ReLU(),
            nn.Linear(256, mutation_dim)
        )
        
        self.mirna_decoder = nn.Sequential(
            nn.Linear(2, 64),
            nn.ReLU(),
            nn.Linear(64, 256),
            nn.ReLU(),
            nn.Linear(256, mirna_dim)
        )
        
        self.embedding_decoder = nn.Sequential(
            nn.Linear(2, 64),
            nn.ReLU(),
            nn.Linear(64, 256),
            nn.ReLU(),
            nn.Linear(256, embedding_dim)
        )
        
    def forward(self, x_mut, x_mirna, x_emb):
        # Encode each modality
        mut_encoded = self.mutation_encoder(x_mut)
        mirna_encoded = self.mirna_encoder(x_mirna)
        emb_encoded = self.embedding_encoder(x_emb)
        
        # Concatenate and fuse
        combined = torch.cat([mut_encoded, mirna_encoded, emb_encoded], dim=1)
        fused = self.fusion(combined)
        
        # Decode each modality
        mut_decoded = self.mutation_decoder(fused)
        mirna_decoded = self.mirna_decoder(fused)
        emb_decoded = self.embedding_decoder(fused)
        
        return fused, (mut_decoded, mirna_decoded, emb_decoded)

def process_and_align_data(mutation_file, mirna_file, embedding_file, id_mapping_file, labels_file):
    """Load and align multi-omics data using patient IDs"""
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
    
    # Load mutation data
    mutations = pd.read_csv(mutation_file)
    mutations = mutations.drop_duplicates(subset='patient_barcode', keep='first')
    mutations = mutations[mutations['patient_barcode'].isin(valid_patients)]
    mutation_features = mutations.set_index('patient_barcode')
    print(f"Loaded mutation data: {mutation_features.shape}")
    
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
    
    # Find common patients
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
    
    # Sort all dataframes by patient ID
    common_patients = sorted(list(common_patients))
    mutation_features = mutation_features.reindex(index=common_patients)
    mirna = mirna.reindex(index=common_patients)
    embeddings_df = embeddings_df.reindex(index=common_patients)
    labels = labels.set_index('patient_barcode').reindex(index=common_patients)
    
    return mutation_features, mirna, embeddings_df, labels

def train_autoencoder(model, train_loader, device, epochs=100):
    """Train autoencoder model"""
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    model.to(device)
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        
        for batch in train_loader:
            if isinstance(model, EarlyFusionAE):
                # Early fusion: concatenate all features
                x_mut, x_mirna, x_emb = [x.to(device) for x in batch]
                x_combined = torch.cat([x_mut, x_mirna, x_emb], dim=1)
                
                encoded, decoded = model(x_combined)
                loss = F.mse_loss(decoded, x_combined)
            else:
                # Late fusion: separate modalities
                x_mut, x_mirna, x_emb = [x.to(device) for x in batch]
                encoded, (mut_decoded, mirna_decoded, emb_decoded) = model(x_mut, x_mirna, x_emb)
                
                loss = (F.mse_loss(mut_decoded, x_mut) + 
                       F.mse_loss(mirna_decoded, x_mirna) + 
                       F.mse_loss(emb_decoded, x_emb))
            
            optimizer.zero_grad()
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
            if isinstance(model, EarlyFusionAE):
                x_mut, x_mirna, x_emb = [x.to(device) for x in batch]
                x_combined = torch.cat([x_mut, x_mirna, x_emb], dim=1)
                encoded, _ = model(x_combined)
            else:
                x_mut, x_mirna, x_emb = [x.to(device) for x in batch]
                encoded, _ = model(x_mut, x_mirna, x_emb)
            
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

def plot_latent_space(latent_features, cluster_labels, high_risk_cluster, title, output_path):
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
    plt.title(title)
    plt.xlabel('Latent Dimension 1')
    plt.ylabel('Latent Dimension 2')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_km_curves(durations, events, cluster_labels, high_risk_cluster, output_path):
    """Create Kaplan-Meier survival plots for risk groups"""
    plt.figure(figsize=(10, 6))
    kmf = KaplanMeierFitter()
    
    # Colors for risk groups
    colors = {'Low-risk': '#66CCCC', 'High-risk': '#FF9999'}
    
    # Plot each risk group
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
    labels_file='/Users/stanleychen/git/Melanoma/data/patient_survival_analysis_filtered_keep_negs.csv'
    
    # Load and process data
    mutation_data, mirna_data, embedding_data, labels = process_and_align_data(
        mutation_file, mirna_file, embedding_file, id_mapping_file, labels_file
    )
    
    # Scale the data
    scaler = StandardScaler()
    mutation_scaled = scaler.fit_transform(mutation_data)
    mirna_scaled = scaler.fit_transform(mirna_data)
    embedding_scaled = scaler.fit_transform(embedding_data)
    
    # Create dataset and dataloader
    dataset = MultiOmicsDataset(mutation_scaled, mirna_scaled, embedding_scaled)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    # Define input dimensions
    input_dims = [mutation_scaled.shape[1], 
                 mirna_scaled.shape[1], 
                 embedding_scaled.shape[1]]
    
    # Device configuration
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"\nUsing device: {device}")
    
    # Train Early Fusion AE
    print("\nTraining Early Fusion Autoencoder...")
    early_fusion_ae = EarlyFusionAE(input_dims)
    early_fusion_ae = train_autoencoder(early_fusion_ae, dataloader, device)
    
    # Train Late Fusion AE
    print("\nTraining Late Fusion Autoencoder...")
    late_fusion_ae = LateFusionAE(input_dims)
    late_fusion_ae = train_autoencoder(late_fusion_ae, dataloader, device)
    
    # Get latent representations
    early_fusion_latent = get_latent_space(early_fusion_ae, dataloader, device)
    late_fusion_latent = get_latent_space(late_fusion_ae, dataloader, device)
    
    # Perform clustering on both latent spaces
    print("\nPerforming clustering...")
    
    # Early fusion clustering
    early_kmeans = KMeans(n_clusters=2, random_state=42)
    early_cluster_labels = early_kmeans.fit_predict(early_fusion_latent)
    
    # Late fusion clustering
    late_kmeans = KMeans(n_clusters=2, random_state=42)
    late_cluster_labels = late_kmeans.fit_predict(late_fusion_latent)
    
    # Calculate survival metrics for both models
    print("\nEarly Fusion Results:")
    early_c_index, early_p_value, early_high_risk = calculate_survival_metrics(
        labels['survival_time_days'].values,
        labels['event'].values,
        early_cluster_labels,
        early_fusion_latent
    )
    
    print("\nLate Fusion Results:")
    late_c_index, late_p_value, late_high_risk = calculate_survival_metrics(
        labels['survival_time_days'].values,
        labels['event'].values,
        late_cluster_labels,
        late_fusion_latent
    )
    
    # Plot results
    plot_latent_space(
        early_fusion_latent,
        early_cluster_labels,
        early_high_risk,
        'Early Fusion AE Clustering',
        'early_fusion_clusters.png'
    )
    
    plot_latent_space(
        late_fusion_latent,
        late_cluster_labels,
        late_high_risk,
        'Late Fusion AE Clustering',
        'late_fusion_clusters.png'
    )
    
    # Plot survival curves
    plot_km_curves(
        labels['survival_time_days'].values,
        labels['event'].values,
        early_cluster_labels,
        early_high_risk,
        'early_fusion_survival.png'
    )
    
    plot_km_curves(
        labels['survival_time_days'].values,
        labels['event'].values,
        late_cluster_labels,
        late_high_risk,
        'late_fusion_survival.png'
    )
    
    # Save cluster assignments
    early_fusion_df = pd.DataFrame({
        'patient_id': labels.index,
        'cluster': early_cluster_labels,
        'risk_group': ['High-risk' if c == early_high_risk else 'Low-risk' 
                      for c in early_cluster_labels],
        'survival_time': labels['survival_time_days'],
        'event': labels['event']
    })
    early_fusion_df.to_csv('early_fusion_clusters.csv')
    
    late_fusion_df = pd.DataFrame({
        'patient_id': labels.index,
        'cluster': late_cluster_labels,
        'risk_group': ['High-risk' if c == late_high_risk else 'Low-risk' 
                      for c in late_cluster_labels],
        'survival_time': labels['survival_time_days'],
        'event': labels['event']
    })
    late_fusion_df.to_csv('late_fusion_clusters.csv')
    
    return {
        'early_fusion': {
            'c_index': early_c_index,
            'p_value': early_p_value,
            'cluster_labels': early_cluster_labels,
            'latent_features': early_fusion_latent
        },
        'late_fusion': {
            'c_index': late_c_index,
            'p_value': late_p_value,
            'cluster_labels': late_cluster_labels,
            'latent_features': late_fusion_latent
        }
    }

if __name__ == "__main__":
    results = main()