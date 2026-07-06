import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import kneighbors_graph

# ---------------------------
# User Inputs: Define File Paths
# ---------------------------
expression_file = '/Users/stanleychen/git/Melanoma/final_data/combined_RNA-seq_data.csv'  # Replace with your actual file path
labels_file = '/Users/stanleychen/git/Melanoma/final_data/final_label.csv'  # Replace with your actual file path
output_file = 'processed_data.pkl'

# ---------------------------
# Step 1: Load Data
# ---------------------------
# Load RNA-Seq expression data
expression_data = pd.read_csv(expression_file, index_col=0)  # Rows: genes, Columns: samples

# Load labels data
labels = pd.read_csv(labels_file)

# Clean column names and identifiers
expression_data.columns = expression_data.columns.str.strip()
labels['patient_barcode'] = labels['patient_barcode'].str.strip()

# ---------------------------
# Step 2: Filter Valid Samples
# ---------------------------
# Remove samples with missing survival_time_days or tumor_stage
labels = labels.dropna(subset=['survival_time_days', 'tumor_stage'])

# Match valid sample IDs between expression and labels
valid_samples = list(set(labels['patient_barcode']).intersection(set(expression_data.columns)))
labels = labels[labels['patient_barcode'].isin(valid_samples)]
expression_data = expression_data[valid_samples]

# ---------------------------
# Step 3: Generate Global Adjacency Matrix
# ---------------------------
def generate_adjacency_matrix(expression_data, method='correlation', k=5):
    if method == 'correlation':
        corr_matrix = expression_data.T.corr()
        adjacency_matrix = (corr_matrix > 0.5).astype(int)
        np.fill_diagonal(adjacency_matrix.values, 0)
        return adjacency_matrix.values
    elif method == 'knn':
        adjacency_matrix = kneighbors_graph(expression_data.T, n_neighbors=k, mode='connectivity').toarray()
        return adjacency_matrix
    else:
        raise ValueError("Invalid method. Use 'correlation' or 'knn'.")

global_graph = generate_adjacency_matrix(expression_data.T, method='correlation')

# ---------------------------
# Step 4: Standardize Tumor Grades
# ---------------------------
def clean_tumor_grade(grade):
    """
    Removes the last character if it is a letter.
    Example: 'Stage IIA' -> 'Stage II', 'Stage III' -> 'Stage III'
    """
    if grade[-1].isalpha():
        return grade[:-1].strip()
    return grade.strip()

# Apply tumor grade cleaning to the labels
labels['tumor_stage'] = labels['tumor_stage'].apply(clean_tumor_grade)

# ---------------------------
# Step 5: Balanced Train, Validation, Test Split
# ---------------------------
# Shuffle data to ensure randomness
labels = labels.sample(frac=1, random_state=42).reset_index(drop=True)

# Perform the split (70% Train, 15% Validation, 15% Test)
train_labels, temp_labels = train_test_split(labels, test_size=0.3, random_state=42)
validation_labels, test_labels = train_test_split(temp_labels, test_size=0.5, random_state=42)

# Function to prepare data splits
def prepare_split(split_labels, expression_data, global_graph):
    split_samples = split_labels['patient_barcode'].tolist()
    split_indices = [valid_samples.index(sample) for sample in split_samples]
    
    split_expression_data = expression_data[split_samples].T
    split_graph = global_graph[np.ix_(split_indices, split_indices)]
    
    return {
        'x_omic_10000': split_expression_data.values,  # Gene expression data
        'x_grph': split_graph,  # Adjacency matrix for this split
        'e': (split_labels['vital_status'] == 'Dead').astype(int).tolist(),  # Event data
        't': split_labels['survival_time_days'].tolist(),  # Survival time
        'g': split_labels['tumor_stage'].tolist()  # Cleaned Tumor grade
    }

# Prepare splits
data_splits = {
    'train': prepare_split(train_labels, expression_data, global_graph),
    'validation': prepare_split(validation_labels, expression_data, global_graph),
    'test': prepare_split(test_labels, expression_data, global_graph)
}

# Save sample splits
cv_splits = {
    'train': train_labels['patient_barcode'].tolist(),
    'validation': validation_labels['patient_barcode'].tolist(),
    'test': test_labels['patient_barcode'].tolist()
}

# ---------------------------
# Step 6: Save to Pickle
# ---------------------------
processed_data = {
    'cv_splits': cv_splits,
    'graph': global_graph,
    'train': data_splits['train'],
    'validation': data_splits['validation'],
    'test': data_splits['test']
}

with open(output_file, 'wb') as f:
    pickle.dump(processed_data, f)

print(f"✅ Processed data saved successfully to {output_file}")
