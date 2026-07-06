import torch
from torch.utils.data import DataLoader
import numpy as np
from sklearn.metrics import (
    accuracy_score, precision_recall_fscore_support,
    confusion_matrix, roc_auc_score
)
from typing import Dict, List, Tuple
import logging
import json
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


class ModelValidator:
    def __init__(
        self,
        model,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu'
    ):
        self.model = model.to(device)
        self.device = device
        self.logger = logging.getLogger(__name__)

    def validate_fold(
        self,
        val_loader: DataLoader,
        fold_num: int
    ) -> Dict[str, float]:
        """Validate a single fold."""
        self.model.eval()
        predictions = []
        true_labels = []
        survival_preds = []
        survival_true = []
        class_probs = []
        
        with torch.no_grad():
            for inputs, survival, labels in val_loader:
                # Move data to device
                inputs = {k: v.to(self.device) for k, v in inputs.items()}
                survival = survival.to(self.device)
                labels = labels.to(self.device).squeeze()
                
                # Forward pass
                survival_pred, class_pred = self.model(inputs)
                
                # Store predictions and true values
                predictions.extend(class_pred.argmax(dim=1).cpu().numpy())
                class_probs.extend(torch.softmax(class_pred, dim=1).cpu().numpy())
                true_labels.extend(labels.cpu().numpy())
                survival_preds.extend(survival_pred.cpu().numpy())
                survival_true.extend(survival.cpu().numpy())
        
        # Calculate metrics
        metrics = self._calculate_metrics(
            np.array(true_labels),
            np.array(predictions),
            np.array(survival_true),
            np.array(survival_preds),
            np.array(class_probs)
        )
        
        return metrics

    def _calculate_metrics(
        self,
        true_labels: np.ndarray,
        predictions: np.ndarray,
        survival_true: np.ndarray,
        survival_preds: np.ndarray,
        class_probs: np.ndarray
    ) -> Dict[str, float]:
        """Calculate comprehensive metrics for both tasks."""
        # Classification metrics
        accuracy = accuracy_score(true_labels, predictions)
        precision, recall, f1, _ = precision_recall_fscore_support(
            true_labels, predictions, average='weighted'
        )
        conf_matrix = confusion_matrix(true_labels, predictions)
        
        # ROC AUC calculation
        try:
            roc_auc = roc_auc_score(true_labels, class_probs, multi_class='ovr')
        except Exception as e:
            self.logger.warning(f"Could not calculate ROC AUC: {str(e)}")
            roc_auc = None
        
        # Survival metrics
        mse = np.mean((survival_true - survival_preds) ** 2)
        mae = np.mean(np.abs(survival_true - survival_preds))
        
        # Correlation between survival predictions and true values
        correlation = np.corrcoef(survival_true.flatten(), survival_preds.flatten())[0, 1]
        
        return {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'roc_auc': roc_auc,
            'survival_mse': mse,
            'survival_mae': mae,
            'survival_correlation': correlation,
            'confusion_matrix': conf_matrix.tolist()
        }

    def _average_metrics(self, fold_metrics: List[Dict]) -> Dict[str, float]:
        """Calculate average metrics across folds."""
        avg_metrics = {}
        metrics_to_avg = [
            'accuracy', 'precision', 'recall', 'f1_score', 
            'survival_mse', 'survival_mae', 'survival_correlation'
        ]
        
        for metric in metrics_to_avg:
            values = [fold[metric] for fold in fold_metrics if metric in fold]
            if values:
                avg_metrics[f'avg_{metric}'] = np.mean(values)
                avg_metrics[f'std_{metric}'] = np.std(values)
        
        return avg_metrics

    def _generate_plots(
        self,
        fold_metrics: List[Dict],
        output_dir: Path
    ) -> None:
        """Generate and save visualization plots."""
        # 1. Confusion Matrix
        plt.figure(figsize=(10, 8))
        avg_conf_matrix = np.mean([np.array(m['confusion_matrix']) 
                                 for m in fold_metrics], axis=0)
        sns.heatmap(
            avg_conf_matrix,
            annot=True,
            fmt='.2f',
            cmap='Blues',
            xticklabels=['I', 'II', 'III', 'IV'],
            yticklabels=['I', 'II', 'III', 'IV']
        )
        plt.title('Average Confusion Matrix Across Folds')
        plt.xlabel('Predicted Stage')
        plt.ylabel('True Stage')
        plt.tight_layout()
        plt.savefig(output_dir / 'confusion_matrix.png')
        plt.close()

        # 2. Metrics Distribution
        for metric in ['accuracy', 'f1_score', 'survival_mae']:
            values = [m[metric] for m in fold_metrics]
            plt.figure(figsize=(8, 6))
            sns.boxplot(data=[values])
            plt.title(f'{metric.replace("_", " ").title()} Distribution Across Folds')
            plt.tight_layout()
            plt.savefig(output_dir / f'{metric}_distribution.png')
            plt.close()

    def validate_model(
        self,
        val_loader: DataLoader,
        model_dir: str,
        output_dir: str
    ) -> Dict[str, float]:
        """Validate model across all folds and save results."""
        model_dir = Path(model_dir)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        fold_metrics = []
        all_predictions = []
        all_true_labels = []
        all_survival_preds = []
        all_survival_true = []
        
        # Load and validate each fold
        for fold_path in sorted(model_dir.glob('fold_*_best_model.pt')):
            fold_num = int(fold_path.stem.split('_')[1])
            self.logger.info(f"Validating fold {fold_num}")
            
            # Load model weights
            checkpoint = torch.load(fold_path)
            self.model.load_state_dict(checkpoint['model_state_dict'])
            
            # Validate fold
            metrics = self.validate_fold(val_loader, fold_num)
            fold_metrics.append(metrics)
            
            # Save fold metrics
            with open(output_dir / f'fold_{fold_num}_metrics.json', 'w') as f:
                json.dump(metrics, f, indent=4)
        
        # Calculate and save average metrics
        avg_metrics = self._average_metrics(fold_metrics)
        with open(output_dir / 'average_metrics.json', 'w') as f:
            json.dump(avg_metrics, f, indent=4)
        
        # Generate and save plots
        self._generate_plots(fold_metrics, output_dir)
        
        # Log summary
        self.logger.info("Validation Summary:")
        for metric, value in avg_metrics.items():
            self.logger.info(f"{metric}: {value:.4f}")
        
        return avg_metrics

    def predict(
        self,
        test_loader: DataLoader,
        model_path: str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Make predictions using the best model."""
        # Load best model
        checkpoint = torch.load(model_path)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.eval()
        
        predictions = []
        survival_preds = []
        
        with torch.no_grad():
            for inputs, *_ in test_loader:
                inputs = {k: v.to(self.device) for k, v in inputs.items()}
                survival_pred, class_pred = self.model(inputs)
                
                predictions.extend(class_pred.argmax(dim=1).cpu().numpy())
                survival_preds.extend(survival_pred.cpu().numpy())
        
        return np.array(predictions), np.array(survival_preds)