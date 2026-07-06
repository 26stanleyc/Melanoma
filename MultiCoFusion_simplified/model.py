import torch
import torch.nn as nn
from typing import Dict, Optional, Tuple

from model_components import (
    FeatureEncoder, CrossModalAttention, 
    ModalityFusion, TaskSpecificHead
)


class FlexibleOmicsModel(nn.Module):
    def __init__(
        self,
        modality_dims: Dict[str, int],
        fusion_dim: int = 256,
        num_classes: int = 4,
        dropout: float = 0.2
    ):
        super().__init__()
        self.modality_dims = modality_dims
        
        # Feature encoders for each modality
        self.encoders = nn.ModuleDict({
            name: FeatureEncoder(dim, dropout=dropout)
            for name, dim in modality_dims.items()
        })
        
        # Cross-modal attention
        self.attention = CrossModalAttention(modality_dims)
        
        # Fusion module
        self.fusion = ModalityFusion(modality_dims, fusion_dim)
        
        # Task heads
        self.task_heads = TaskSpecificHead(fusion_dim, num_classes)
        
    def forward(
        self, 
        inputs: Dict[str, Optional[torch.Tensor]]
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        # Filter out None inputs and their corresponding modalities
        valid_inputs = {
            name: tensor 
            for name, tensor in inputs.items() 
            if tensor is not None and name in self.encoders
        }
        
        if not valid_inputs:
            raise ValueError("No valid inputs provided")
        
        # Encode each modality
        encoded = {
            name: self.encoders[name](tensor)
            for name, tensor in valid_inputs.items()
        }
        
        # Cross-modal attention
        attended = self.attention(encoded)
        
        # Fuse modalities
        fused = self.fusion(attended)
        
        # Get task-specific predictions
        survival_pred, class_pred = self.task_heads(fused)
        
        return survival_pred, class_pred