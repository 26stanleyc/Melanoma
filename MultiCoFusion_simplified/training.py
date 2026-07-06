import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np
from sklearn.model_selection import KFold
from typing import Dict, List, Tuple
import logging
from tqdm import tqdm
from pathlib import Path
import json


class Trainer:
    def __init__(
        self,
        model: nn.Module,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu',
        lr: float = 1e-5,  # Reduced learning rate
        weight_decay: float = 1e-4,  # Increased weight decay
        alpha: float = 0.5
    ):
        self.model = model.to(device)
        self.device = device
        self.optimizer = torch.optim.Adam(
            model.parameters(),
            lr=lr,
            weight_decay=weight_decay
        )
        self.alpha = alpha
        
        # Initialize logger
        self.logger = logging.getLogger(__name__)
        
        # Calculate class weights based on distribution
        class_counts = np.array([277, 12, 0, 21])  # Your class distribution
        total_samples = class_counts.sum()
        # More stable class weight calculation
        class_weights = total_samples / (class_counts + 1e-6) / len(class_counts)
        class_weights = torch.FloatTensor(class_weights).to(device)
        
        # Loss functions with stability modifications
        self.survival_criterion = nn.MSELoss(reduction='mean')
        self.classification_criterion = nn.CrossEntropyLoss(
            weight=class_weights,
            reduction='mean'
        )
        
        # More conservative scheduler
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer,
            mode='min',
            factor=0.2,
            patience=10,
            min_lr=1e-7,
            verbose=True
        )

    def train_step(
        self,
        inputs: Dict[str, torch.Tensor],
        survival_time: torch.Tensor,
        stage_labels: torch.Tensor
    ) -> Dict[str, float]:
        self.model.train()
        
        # Move data to device and handle NaN values
        inputs = {k: torch.nan_to_num(v, nan=0.0).to(self.device) 
                 for k, v in inputs.items()}
        survival_time = torch.nan_to_num(survival_time, nan=0.0).to(self.device)
        stage_labels = stage_labels.to(self.device).squeeze()
        
        # Forward pass
        survival_pred, class_pred = self.model(inputs)
        
        # Add small epsilon to avoid log(0)
        eps = 1e-7
        survival_pred = torch.clamp(survival_pred, min=-1e6, max=1e6)
        
        # Calculate losses with stability checks
        try:
            survival_loss = self.survival_criterion(survival_pred, survival_time)
            classification_loss = self.classification_criterion(class_pred, stage_labels)
            
            # Check for NaN losses
            if torch.isnan(survival_loss):
                self.logger.warning("NaN detected in survival loss")
                survival_loss = torch.tensor(0.0, device=self.device)
            if torch.isnan(classification_loss):
                self.logger.warning("NaN detected in classification loss")
                classification_loss = torch.tensor(0.0, device=self.device)
            
            # Calculate dynamic alpha for loss weighting
            current_epoch = getattr(self, 'current_epoch', 0)
            warmup_epochs = 10
            dynamic_alpha = min(self.alpha, self.alpha * (current_epoch / warmup_epochs))
            
            # Scale survival loss based on its magnitude compared to classification loss
            loss_scale = torch.log1p(classification_loss) / torch.log1p(survival_loss + eps)
            scaled_survival_loss = survival_loss * loss_scale.detach()
            
            total_loss = dynamic_alpha * scaled_survival_loss + (1 - dynamic_alpha) * classification_loss
            
            # Backward pass with gradient clipping
            self.optimizer.zero_grad()
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()
            
            # Calculate metrics
            with torch.no_grad():
                pred = class_pred.argmax(dim=1)
                accuracy = (pred == stage_labels).float().mean().item()
            
            metrics = {
                'total_loss': total_loss.item(),
                'survival_loss': survival_loss.item(),
                'scaled_survival_loss': scaled_survival_loss.item(),
                'classification_loss': classification_loss.item(),
                'accuracy': accuracy,
                'dynamic_alpha': dynamic_alpha,
                'loss_scale': loss_scale.item()
            }
            
            # Check for NaN metrics
            metrics = {k: v if not np.isnan(v) else 0.0 for k, v in metrics.items()}
            
            return metrics
            
        except RuntimeError as e:
            self.logger.error(f"Runtime error in training step: {str(e)}")
            return {
                'total_loss': 0.0,
                'survival_loss': 0.0,
                'classification_loss': 0.0,
                'accuracy': 0.0,
                'dynamic_alpha': 0.0,
                'loss_scale': 0.0
            }
        
    @torch.no_grad()
    def validate(self, val_loader: DataLoader) -> Dict[str, float]:
        self.model.eval()
        total_class_loss = 0.0
        total_survival_loss = 0.0
        total_correct = 0
        total_samples = 0
        survival_preds = []
        survival_true = []
        
        with torch.no_grad():
            for batch in val_loader:
                inputs, survival_time, stage_labels = batch
                
                # Move to device
                inputs = {k: v.to(self.device) for k, v in inputs.items()}
                survival_time = survival_time.to(self.device)
                stage_labels = stage_labels.to(self.device).squeeze()
                
                # Forward pass
                survival_pred, class_pred = self.model(inputs)
                
                # Calculate losses
                survival_loss = self.survival_criterion(survival_pred, survival_time)
                classification_loss = self.classification_criterion(class_pred, stage_labels)
                
                # Accumulate metrics
                total_class_loss += classification_loss.item()
                total_survival_loss += survival_loss.item()
                pred = class_pred.argmax(dim=1)
                total_correct += (pred == stage_labels).sum().item()
                total_samples += stage_labels.size(0)
                
                # Store survival predictions for correlation calculation
                survival_preds.extend(survival_pred.cpu().numpy())
                survival_true.extend(survival_time.cpu().numpy())
        
        # Calculate correlation for survival predictions
        survival_correlation = np.corrcoef(
            np.array(survival_preds).flatten(),
            np.array(survival_true).flatten()
        )[0, 1]
        
        metrics = {
            'val_class_loss': total_class_loss / len(val_loader),
            'val_survival_loss': total_survival_loss / len(val_loader),
            'val_accuracy': total_correct / total_samples,
            'val_survival_corr': survival_correlation
        }
        
        return metrics

    def _init_weights(self, module):
        """Initialize weights with stability in mind."""
        if isinstance(module, (nn.Linear, nn.Embedding)):
            # Use smaller initialization values
            nn.init.xavier_uniform_(module.weight, gain=0.5)
            if module.bias is not None:
                nn.init.constant_(module.bias, 0.0)
    def get_balanced_sampler(self, dataset_subset):
        # Get labels from the underlying dataset
        if isinstance(dataset_subset, torch.utils.data.Subset):
            # Convert labels to numpy array and then index
            full_labels = np.array(dataset_subset.dataset.stage_labels)
            labels = full_labels[dataset_subset.indices]
        else:
            # If it's the full dataset
            labels = np.array(dataset_subset.stage_labels)

        class_counts = np.bincount(labels)
        weights = 1. / class_counts
        sample_weights = weights[labels]
        sample_weights = torch.from_numpy(sample_weights).float()
        sampler = torch.utils.data.WeightedRandomSampler(
            weights=sample_weights,
            num_samples=len(sample_weights),
            replacement=True
        )
        return sampler
    
    def train_with_cross_validation(
        self,
        dataset,
        save_dir: str,
        n_splits: int = 5,
        num_epochs: int = 100,
        batch_size: int = 32,
        early_stopping_patience: int = 15
    ) -> List[Dict[str, float]]:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        
        kfold = KFold(n_splits=n_splits, shuffle=True)
        cv_results = []
        
        for fold, (train_idx, val_idx) in enumerate(kfold.split(dataset)):
            self.logger.info(f"Starting fold {fold+1}/{n_splits}")
            
            # Create data loaders for this fold
            train_subsampler = torch.utils.data.SubsetRandomSampler(train_idx)
            val_subsampler = torch.utils.data.SubsetRandomSampler(val_idx)
            train_subset = torch.utils.data.Subset(dataset, train_idx)

            train_sampler = self.get_balanced_sampler(train_subset)
            train_loader = DataLoader(
                dataset,
                batch_size=batch_size,
                sampler=train_sampler,
                num_workers=4,
                pin_memory=True
            )
            val_loader = DataLoader(
                dataset,
                batch_size=batch_size,
                sampler=val_subsampler,
                num_workers=4,
                pin_memory=True
            )
            
            # Training loop for this fold
            best_val_loss = float('inf')
            patience_counter = 0
            fold_results = []
            
            # Initialize model weights for this fold
            self._init_weights(self.model)
            
            for epoch in tqdm(range(num_epochs), desc=f"Fold {fold+1}"):
                self.current_epoch = epoch  # Set current epoch for dynamic alpha
                
                # Training phase
                train_metrics = []
                for batch in train_loader:
                    metrics = self.train_step(*batch)
                    train_metrics.append(metrics)
                
                # Calculate average training metrics
                avg_train_metrics = {
                    key: np.mean([m[key] for m in train_metrics])
                    for key in train_metrics[0]
                }
                
                # Validation phase
                val_metrics = self.validate(val_loader)
                
                # Calculate combined validation loss for scheduler and early stopping
                val_total_loss = val_metrics['val_class_loss'] + val_metrics['val_survival_loss']
                
                # Learning rate scheduling
                self.scheduler.step(val_total_loss)
                
                # Save metrics
                epoch_metrics = {
                    'epoch': epoch + 1,
                    **avg_train_metrics,
                    **val_metrics
                }
                fold_results.append(epoch_metrics)
                
                # Early stopping check
                if val_total_loss < best_val_loss:
                    best_val_loss = val_total_loss
                    patience_counter = 0
                    
                    # Save best model for this fold
                    torch.save({
                        'model_state_dict': self.model.state_dict(),
                        'optimizer_state_dict': self.optimizer.state_dict(),
                        'val_metrics': val_metrics,
                        'epoch': epoch
                    }, save_dir / f'fold_{fold+1}_best_model.pt')
                else:
                    patience_counter += 1
                
                # Log progress
                self.logger.info(
                    f"Epoch {epoch+1}/{num_epochs} - "
                    f"Train Loss: {avg_train_metrics['total_loss']:.4f}, "
                    f"Val Class Loss: {val_metrics['val_class_loss']:.4f}, "
                    f"Val Survival Loss: {val_metrics['val_survival_loss']:.4f}, "
                    f"Val Acc: {val_metrics['val_accuracy']:.4f}, "
                    f"Val Survival Corr: {val_metrics['val_survival_corr']:.4f}"
                )
                
                if patience_counter >= early_stopping_patience:
                    self.logger.info(f"Early stopping triggered at epoch {epoch+1}")
                    break
            
            # Save fold results
            cv_results.append({
                'fold': fold + 1,
                'best_val_loss': best_val_loss,
                'results': fold_results
            })
            
            # Save fold results to file
            with open(save_dir / f'fold_{fold+1}_results.json', 'w') as f:
                json.dump(fold_results, f, indent=4)
        
        # Save overall CV results
        with open(save_dir / 'cv_results.json', 'w') as f:
            json.dump(cv_results, f, indent=4)
        
        return cv_results