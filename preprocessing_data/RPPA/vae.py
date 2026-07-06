import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import os
from datetime import datetime
from torch.utils.data import Dataset, DataLoader

class RPPADataset(Dataset):
    def __init__(self, data_path):
        # Read and process RPPA data
        self.data = pd.read_csv(data_path, index_col=0)
        self.data = self.data.astype(float)
        # Normalize the data
        self.data = (self.data - self.data.mean()) / self.data.std()
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        return torch.FloatTensor(self.data.iloc[idx].values)

class RPPAVAE(nn.Module):
    def __init__(self, input_dim, latent_dim=32, hidden_dim=64):
        super(RPPAVAE, self).__init__()
        
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU()
        )
        
        # Latent space
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_var = nn.Linear(hidden_dim, latent_dim)
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim)
        )
        
    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_var(h)
    
    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode(self, z):
        return self.decoder(z)
    
    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        return self.decode(z), mu, log_var

def get_latent_representations(model, data_path):
    dataset = RPPADataset(data_path)
    loader = DataLoader(dataset, batch_size=32, shuffle=False)
    device = next(model.parameters()).device
    
    latent_vectors = []
    model.eval()
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            mu, _ = model.encode(batch)
            latent_vectors.append(mu.cpu())
    
    return torch.cat(latent_vectors, dim=0)

def run_rppa_vae(data_path, output_dir="vae_output", latent_dim=32, hidden_dim=64, epochs=100):
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Initialize dataset and dataloader
    dataset = RPPADataset(data_path)
    input_dim = dataset[0].shape[0]
    train_loader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    # Initialize model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = RPPAVAE(
        input_dim=input_dim,
        latent_dim=latent_dim,
        hidden_dim=hidden_dim
    ).to(device)
    
    # Training
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    
    # Training log
    training_log = []
    
    model.train()
    for epoch in range(epochs):
        total_loss = 0
        for batch in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad()
            
            recon_batch, mu, log_var = model(batch)
            
            # Compute losses
            recon_loss = F.mse_loss(recon_batch, batch, reduction='sum')
            kl_loss = -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())
            loss = recon_loss + kl_loss
            
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
        avg_loss = total_loss / len(train_loader.dataset)
        training_log.append(avg_loss)
        
        if (epoch + 1) % 10 == 0:
            print(f'Epoch [{epoch+1}/{epochs}], Average Loss: {avg_loss:.4f}')
    
    # Save the model
    model_path = os.path.join(output_dir, f"rppa_vae_model_{timestamp}.pt")
    torch.save({
        'model_state_dict': model.state_dict(),
        'latent_dim': latent_dim,
        'hidden_dim': hidden_dim,
        'input_dim': input_dim,
        'training_log': training_log
    }, model_path)
    print(f"Model saved to {model_path}")
    
    # Save latent representations
    latent_reps = get_latent_representations(model, data_path)
    latent_path = os.path.join(output_dir, f"latent_representations_{timestamp}.pt")
    torch.save(latent_reps, latent_path)
    print(f"Latent representations saved to {latent_path}")
    
    # Save training log
    log_path = os.path.join(output_dir, f"training_log_{timestamp}.csv")
    pd.DataFrame({'loss': training_log}).to_csv(log_path)
    print(f"Training log saved to {log_path}")
    
    return model, latent_reps, training_log

def load_rppa_vae(model_path):
    checkpoint = torch.load(model_path)
    model = RPPAVAE(
        input_dim=checkpoint['input_dim'],
        latent_dim=checkpoint['latent_dim'],
        hidden_dim=checkpoint['hidden_dim']
    )
    model.load_state_dict(checkpoint['model_state_dict'])
    return model, checkpoint.get('training_log', None)

if __name__ == "__main__":
    # Example usage
    data_path = "/Users/stanleychen/git/Melanoma/final_data/combined_RPPA_data.csv"
    
    # Train and save model
    model, latent_reps, training_log = run_rppa_vae(
        data_path,
        output_dir="vae_output",
        latent_dim=32,
        hidden_dim=64,
        epochs=100
    )
    
    # Later, to load the model:
    # model_path = "vae_output/rppa_vae_model_TIMESTAMP.pt"
    # loaded_model, training_log = load_rppa_vae(model_path)