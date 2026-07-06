import torch
import torch.nn as nn
import torch.nn.functional as F


class FeatureEncoder(nn.Module):
    """
    Light encoder for already reduced features.
    Each modality keeps its own dimension until fusion.
    """
    def __init__(self, input_dim: int, factor: float = 1.5, dropout: float = 0.2):
        super().__init__()
        # Each modality keeps its native dimension but gets enhanced
        hidden_dim = int(input_dim * factor)  # Allow slight expansion for feature learning
        
        self.encoder = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, input_dim),  # Return to original dimension
            nn.LayerNorm(input_dim)
        )
        
    def forward(self, x):
        return self.encoder(x)


class CrossModalAttention(nn.Module):
    """
    Attention mechanism that can handle different input dimensions.
    Projects each modality to a common space for attention.
    """
    def __init__(self, modal_dims: dict, common_dim: int = 128):
        super().__init__()
        self.common_dim = common_dim
        
        # Projection layers for each modality
        self.projections = nn.ModuleDict({
            name: nn.Linear(dim, common_dim)
            for name, dim in modal_dims.items()
        })
        
        # Multi-head attention in common space
        self.attention = nn.MultiheadAttention(
            common_dim, 
            num_heads=8, 
            batch_first=True
        )
        
        # Output projections back to original spaces
        self.output_projections = nn.ModuleDict({
            name: nn.Linear(common_dim, dim)
            for name, dim in modal_dims.items()
        })
        
    def forward(self, features_dict):
        if len(features_dict) == 1:
            # Single modality case
            modality, features = next(iter(features_dict.items()))
            return {modality: features}
            
        # Project all features to common space
        projected_features = {
            name: self.projections[name](features)
            for name, features in features_dict.items()
        }
        
        # Stack for attention
        stacked = torch.stack(list(projected_features.values()), dim=1)
        
        # Self-attention in common space
        attended, _ = self.attention(stacked, stacked, stacked)
        
        # Project back to original spaces
        outputs = {}
        for idx, (name, _) in enumerate(features_dict.items()):
            modality_attended = attended[:, idx]
            outputs[name] = self.output_projections[name](modality_attended)
            
        return outputs


class ModalityFusion(nn.Module):
    """
    Fuses different modalities while preserving their individual characteristics.
    """
    def __init__(self, modal_dims: dict, fusion_dim: int = 256):
        super().__init__()
        self.fusion_dim = fusion_dim
        
        # Projections to fusion space
        self.projections = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(dim, fusion_dim),
                nn.LayerNorm(fusion_dim),
                nn.ReLU()
            )
            for name, dim in modal_dims.items()
        })
        
        # Gating mechanism
        self.gates = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(dim, 1),
                nn.Sigmoid()
            )
            for name, dim in modal_dims.items()
        })
        
    def forward(self, features_dict):
        if len(features_dict) == 1:
            # Single modality case
            modality, features = next(iter(features_dict.items()))
            return self.projections[modality](features)
        
        # Project each modality
        projected = {
            name: self.projections[name](features)
            for name, features in features_dict.items()
        }
        
        # Calculate gates
        gates = {
            name: self.gates[name](features)
            for name, features in features_dict.items()
        }
        
        # Weighted combination
        combined = torch.zeros(
            features_dict[next(iter(features_dict))].size(0),
            self.fusion_dim
        ).to(next(iter(features_dict.values())).device)
        
        for name in features_dict:
            combined += projected[name] * gates[name]
            
        return combined


class TaskSpecificHead(nn.Module):
    """Output heads for survival and classification tasks."""
    def __init__(self, input_dim: int, num_classes: int):
        super().__init__()
        
        self.shared = nn.Sequential(
            nn.Linear(input_dim, input_dim * 2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(input_dim * 2, input_dim)
        )
        
        self.survival = nn.Sequential(
            nn.Linear(input_dim, input_dim // 2),
            nn.ReLU(),
            nn.Linear(input_dim // 2, 1)
        )
        
        self.classification = nn.Sequential(
            nn.Linear(input_dim, input_dim // 2),
            nn.ReLU(),
            nn.Linear(input_dim // 2, num_classes)
        )
    
    def forward(self, x):
        shared_features = self.shared(x)
        survival = self.survival(shared_features)
        classification = self.classification(shared_features)
        return survival, classification