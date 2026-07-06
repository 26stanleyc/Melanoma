import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
import numpy as np

class RPPAVAE(nn.Module):
    def __init__(self, input_dim, latent_dim=2, hidden_dims=[128, 64], beta=1.0):
        super(RPPAVAE, self).__init__()
        self.latent_dim = latent_dim
        self.beta = beta  # KL weight for beta-VAE functionality

        # Encoder layers
        encoder_layers = []
        prev_dim = input_dim
        for hidden_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.LeakyReLU(),
                nn.Dropout(0.2)
            ])
            prev_dim = hidden_dim
        self.encoder = nn.Sequential(*encoder_layers)

        # Latent space projections
        self.mu = nn.Linear(hidden_dims[-1], latent_dim)
        self.log_var = nn.Linear(hidden_dims[-1], latent_dim)

        # Decoder layers
        decoder_layers = []
        prev_dim = latent_dim
        for hidden_dim in reversed(hidden_dims):
            decoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.LeakyReLU(),
                nn.Dropout(0.2)
            ])
            prev_dim = hidden_dim
        decoder_layers.append(nn.Linear(hidden_dims[0], input_dim))
        self.decoder = nn.Sequential(*decoder_layers)

    def encode(self, x):
        x = self.encoder(x)
        return self.mu(x), self.log_var(x)

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

class RPPAVAETrainer:
    def __init__(self, model, learning_rate=1e-3, device='cuda' if torch.cuda.is_available() else 'cpu'):
        self.model = model.to(device)
        self.device = device
        self.optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        
    def compute_loss(self, x, x_recon, mu, log_var):
        # Reconstruction loss (using mean squared error for RPPA data)
        recon_loss = F.mse_loss(x_recon, x, reduction='sum')
        
        # KL divergence loss
        kl_loss = -0.5 * torch.sum(1 + log_var - mu.pow(2) - log_var.exp())
        
        # Total loss with beta weighting
        total_loss = recon_loss + self.model.beta * kl_loss
        
        return total_loss, recon_loss, kl_loss

    def train_epoch(self, train_loader):
        self.model.train()
        total_loss = 0
        recon_losses = 0
        kl_losses = 0
        
        for batch_idx, data in enumerate(train_loader):
            data = data[0].to(self.device)
            self.optimizer.zero_grad()
            
            # Forward pass
            recon_batch, mu, log_var = self.model(data)
            
            # Compute losses
            loss, recon_loss, kl_loss = self.compute_loss(data, recon_batch, mu, log_var)
            
            # Backward pass and optimization
            loss.backward()
            self.optimizer.step()
            
            total_loss += loss.item()
            recon_losses += recon_loss.item()
            kl_losses += kl_loss.item()
            
        avg_loss = total_loss / len(train_loader.dataset)
        avg_recon = recon_losses / len(train_loader.dataset)
        avg_kl = kl_losses / len(train_loader.dataset)
        
        return avg_loss, avg_recon, avg_kl

def prepare_rppa_data(data, batch_size=32):
    """Prepare RPPA data for VAE training"""
    # Convert to torch tensor and normalize
    data_tensor = torch.FloatTensor(data)
    data_tensor = (data_tensor - data_tensor.mean(dim=0)) / data_tensor.std(dim=0)
    
    # Create dataloader
    dataset = TensorDataset(data_tensor)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    return loader

# Example usage
def train_rppa_vae(rppa_data, n_epochs=100, latent_dim=2, batch_size=32):
    # Prepare data
    train_loader = prepare_rppa_data(rppa_data, batch_size)
    input_dim = rppa_data.shape[1]
    
    # Initialize model and trainer
    model = RPPAVAE(input_dim=input_dim, latent_dim=latent_dim)
    trainer = RPPAVAETrainer(model)
    
    # Training loop
    for epoch in range(n_epochs):
        avg_loss, avg_recon, avg_kl = trainer.train_epoch(train_loader)
        if epoch % 10 == 0:
            print(f'Epoch {epoch}: Loss = {avg_loss:.4f} '
                  f'(Recon = {avg_recon:.4f}, KL = {avg_kl:.4f})')
    
    return model