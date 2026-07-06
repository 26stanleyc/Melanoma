import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from deap import base, creator, tools, algorithms
from lifelines.utils import concordance_index  # Add this import at the top


# Load feature data
features_df = pd.read_csv('/Users/stanleychen/git/Melanoma/final_data/normalized_combined_miRNA.csv', index_col=0)

# Transpose the data so that patient identifiers are rows
features_df = features_df.T

# Reset the index to make patient identifiers a column
features_df = features_df.reset_index()

# Rename the identifier column to match the label file
features_df.rename(columns={'index': 'patient_barcode'}, inplace=True)

# Load survival data
survival_df = pd.read_csv('/Users/stanleychen/git/Melanoma/final_data/filtered_patient_data_all_v2.csv')

# Merge data on 'PatientID'
data = pd.merge(features_df, survival_df, on='patient_barcode')
data['vital_status'] = data['vital_status'].map({'Dead': 1, 'Alive': 0})

# print(data.shape)  # Should return a non-zero shape
# print(data.head())

# Extract features, survival time, and event status
X = data.drop(columns=['patient_barcode', 'survival_time_days', 'vital_status', 'survival_time_source','tumor_stage','grade'])
y_time = data['survival_time_days']
y_event = data['vital_status']

print(X.shape)  # Should return a non-zero shape
print(X.head())

feature_names = X.columns.tolist()

# Standardize Features
scaler = StandardScaler()
X = scaler.fit_transform(X)

# Define DeepSurv Model
class DeepSurv(nn.Module):
    def __init__(self, input_dim):
        super(DeepSurv, self).__init__()
        self.fc1 = nn.Linear(input_dim, 100)
        self.fc2 = nn.Linear(100, 50)
        self.fc3 = nn.Linear(50, 1)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = self.fc3(x)
        return x

# Define Loss Function
def neg_log_partial_likelihood(preds, y_time, y_event):
    hazard_ratio = torch.exp(preds)
    log_risk = torch.log(torch.cumsum(hazard_ratio, dim=0))
    uncensored_likelihood = preds - log_risk
    censored_likelihood = uncensored_likelihood * y_event
    neg_likelihood = -torch.sum(censored_likelihood)
    return neg_likelihood

# Define Fitness Function
def evaluate(individual):
    selected_features = [index for index, value in enumerate(individual) if value == 1]
    if len(selected_features) == 0:
        return 0,

    # Select features
    X_selected = X[:, selected_features]

    # Split data
    X_train, X_val, y_train_time, y_val_time, y_train_event, y_val_event = train_test_split(
        X_selected, y_time, y_event, test_size=0.2, random_state=42
    )

    # Convert to PyTorch tensors
    X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
    y_train_time_tensor = torch.tensor(y_train_time.values, dtype=torch.float32)
    y_train_event_tensor = torch.tensor(y_train_event.values, dtype=torch.float32)
    X_val_tensor = torch.tensor(X_val, dtype=torch.float32)
    y_val_time_tensor = torch.tensor(y_val_time.values, dtype=torch.float32)
    y_val_event_tensor = torch.tensor(y_val_event.values, dtype=torch.float32)

    # Create DataLoader with drop_last=True to ensure consistent batch sizes
    train_dataset = TensorDataset(X_train_tensor, y_train_time_tensor, y_train_event_tensor)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, drop_last=True)

    # Initialize DeepSurv model
    model = DeepSurv(input_dim=X_selected.shape[1])
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    # Train DeepSurv model
    for epoch in range(50):  # 50 epochs
        for X_batch, y_time_batch, y_event_batch in train_loader:
            optimizer.zero_grad()
            
            # Forward pass
            preds = model(X_batch).squeeze()
            
            # Skip if batch is empty
            if preds.nelement() == 0:
                continue
                
            # Make sure preds is 1D
            if len(preds.shape) == 0:
                preds = preds.unsqueeze(0)
            
            # Sort by descending survival time
            order = torch.argsort(y_time_batch, descending=True)
            
            # Apply sorting
            preds = preds[order]
            y_time_batch = y_time_batch[order]
            y_event_batch = y_event_batch[order]

            # Calculate loss
            loss = neg_log_partial_likelihood(preds, y_time_batch, y_event_batch)
            loss.backward()
            optimizer.step()

    # Evaluate on Validation Data
    model.eval()
    with torch.no_grad():
        val_preds = model(X_val_tensor).squeeze()
        
        # Make sure predictions are 1D
        if len(val_preds.shape) == 0:
            val_preds = val_preds.unsqueeze(0)
        
        # Sort by descending survival time
        val_order = torch.argsort(y_val_time_tensor, descending=True)
        val_preds = val_preds[val_order]
        y_val_time_sorted = y_val_time_tensor[val_order]
        y_val_event_sorted = y_val_event_tensor[val_order]

        # Calculate concordance index
        c_index = concordance_index(
            y_val_time_sorted.numpy(),
            -val_preds.numpy(),  # Negative because higher predictions should mean higher risk
            y_val_event_sorted.numpy()
        )

    return c_index,


# Setup Genetic Algorithm
creator.create('FitnessMax', base.Fitness, weights=(1.0,))
creator.create('Individual', list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register('attr_bool', np.random.randint, 2)
toolbox.register('individual', tools.initRepeat, creator.Individual, toolbox.attr_bool, n=X.shape[1])
toolbox.register('population', tools.initRepeat, list, toolbox.individual)
toolbox.register('mate', tools.cxTwoPoint)
toolbox.register('mutate', tools.mutFlipBit, indpb=0.05)
toolbox.register('select', tools.selTournament, tournsize=3)
toolbox.register('evaluate', evaluate)

# Run Genetic Algorithm
population = toolbox.population(n=50)
for gen in range(40):  # 40 generations
    offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.2)
    fits = toolbox.map(toolbox.evaluate, offspring)
    for fit, ind in zip(fits, offspring):
        ind.fitness.values = fit
    population = toolbox.select(offspring, k=len(population))

# Extract Best Features
best_individual = tools.selBest(population, k=1)[0]
selected_features = [index for index, value in enumerate(best_individual) if value == 1]

# Save Selected Features using the stored feature names
selected_feature_names = [feature_names[i] for i in selected_features]
pd.DataFrame(selected_feature_names, columns=['Selected_Features']).to_csv('selected_features.csv', index=False)