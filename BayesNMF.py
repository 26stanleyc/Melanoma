import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from tqdm import tqdm

def bayesnmf(X, k, alpha=0.1, beta=0.1, max_iter=1000, tol=1e-6):
    n, m = X.shape
    W = np.random.gamma(alpha, 1/beta, (n, k))
    H = np.random.gamma(alpha, 1/beta, (k, m))
    
    for _ in range(max_iter):
        W_old, H_old = W.copy(), H.copy()
        
        H = H * (W.T @ X) / (W.T @ W @ H + 1e-10)
        W = W * (X @ H.T) / (W @ H @ H.T + 1e-10)
        
        if np.max(np.abs(W - W_old)) < tol and np.max(np.abs(H - H_old)) < tol:
            break
    
    return W, H

def consensus_clustering(data, max_k, n_iter=100):
    n_samples = data.shape[0]
    consensus_matrix = np.zeros((n_samples, n_samples))
    
    for _ in tqdm(range(n_iter), desc="Consensus Clustering"):
        subsample = np.random.choice(n_samples, size=int(0.8*n_samples), replace=False)
        subdata = data[subsample]
        
        W, H = bayesnmf(subdata.T, max_k)
        labels = np.argmax(W, axis=1)
        
        for i in range(len(labels)):
            for j in range(i+1, len(labels)):
                if labels[i] == labels[j]:
                    consensus_matrix[subsample[i], subsample[j]] += 1
                    consensus_matrix[subsample[j], subsample[i]] += 1
    
    consensus_matrix /= n_iter
    return consensus_matrix

# Load and preprocess data
data = pd.read_csv('final_data/normalized_data/normalized_rna_seq_data.csv', index_col=0)

# Select top 25% highly variable genes
gene_var = data.var()
top_genes = gene_var.nlargest(int(0.25 * len(gene_var))).index
R = data[top_genes]

# Transform data
R_star = R - R.median()

# Compute distance matrix
corr_matrix = R_star.T.corr(method='spearman')
distance_matrix = 1 - corr_matrix

# Consensus Clustering
max_k = 10  # You may need to adjust this
print("Starting Consensus Clustering...")
consensus_matrix = consensus_clustering(R_star.values, max_k)

# Determine optimal number of clusters
cophenetic_corr = []
print("Determining optimal number of clusters...")
for k in tqdm(range(2, max_k+1), desc="KMeans Clustering"):
    kmeans = KMeans(n_clusters=k, n_init=10)
    labels = kmeans.fit_predict(consensus_matrix)
    link = linkage(consensus_matrix, method='average')
    c, _ = spearmanr(consensus_matrix.flatten(), dendrogram(link, no_plot=True)['icoord'])
    cophenetic_corr.append(c)

optimal_k = np.argmax(cophenetic_corr) + 2

# Final clustering
print("Performing final clustering...")
final_labels = KMeans(n_clusters=optimal_k, n_init=10).fit_predict(consensus_matrix)

# Visualize results
plt.figure(figsize=(10, 5))
plt.plot(range(2, max_k+1), cophenetic_corr)
plt.xlabel('Number of clusters (k)')
plt.ylabel('Cophenetic correlation')
plt.title('Optimal number of clusters')
plt.show()

plt.figure(figsize=(10, 10))
plt.imshow(consensus_matrix, cmap='viridis')
plt.title('Consensus Matrix')
plt.colorbar()
plt.show()

print(f"Optimal number of clusters: {optimal_k}")
print("Cluster assignments:", final_labels)