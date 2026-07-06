import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from lifelines.utils import concordance_index
from tqdm import tqdm
from networks import RNASeqNet, GSNN, SGCN

class Trainer:
    def __init__(self, model: nn.Module, device='cuda' if torch.cuda.is_available() else 'cpu', lr=1e-4, weight_decay=1e-4, alpha=0.5):
        self.model = model.to(device)
        self.device = device
        self.optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
        self.alpha = alpha
        
        self.gsnn = GSNN(input_dim=10000, omic_dim=128)
        self.sgcn = SGCN(hidden_dim=128, output_dim=64)
        
        self.survival_head = nn.Linear(64, 1)
        self.classification_head = nn.Linear(64, 3)
        
        self.survival_criterion = nn.MSELoss(reduction='mean')
        self.classification_criterion = nn.CrossEntropyLoss()
        self.scheduler = optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, mode='min', factor=0.2, patience=10, min_lr=1e-7, verbose=True)
    
    def train_step(self, batch: dict):
        self.model.train()
        self.optimizer.zero_grad()
        
        x_omic = batch['x_omic'].to(self.device)
        time = batch['time'].to(self.device)
        event = batch['event'].to(self.device)
        grade = batch['grade'].to(self.device)
        
        # GSNN for Feature Extraction
        features = self.gsnn(x_omic)
        
        # SGCN for Graph-based Feature Interaction
        predictions = self.sgcn(features)
        
        # Two Heads for Prediction
        survival_pred = self.survival_head(predictions)
        classification_pred = self.classification_head(predictions)
        
        survival_loss = self.survival_criterion(survival_pred.squeeze(), time)
        classification_loss = self.classification_criterion(classification_pred, grade)
        total_loss = self.alpha * survival_loss + (1 - self.alpha) * classification_loss
        
        total_loss.backward()
        nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
        self.optimizer.step()
        
        acc = (classification_pred.argmax(dim=1) == grade).float().mean().item()
        return {
            'total_loss': total_loss.item(),
            'survival_loss': survival_loss.item(),
            'classification_loss': classification_loss.item(),
            'accuracy': acc
        }
    
    @torch.no_grad()
    def validate(self, dataloader: DataLoader):
        self.model.eval()
        survival_losses = []
        classification_losses = []
        accuracies = []
        survival_preds = []
        survival_true = []
        
        for batch in dataloader:
            x_omic = batch['x_omic'].to(self.device)
            time = batch['time'].to(self.device)
            event = batch['event'].to(self.device)
            grade = batch['grade'].to(self.device)
            
            # GSNN for Feature Extraction
            features = self.gsnn(x_omic)
            
            # SGCN for Graph-based Feature Interaction
            predictions = self.sgcn(features)
            
            # Two Heads for Prediction
            survival_pred = self.survival_head(predictions)
            classification_pred = self.classification_head(predictions)
            
            survival_loss = self.survival_criterion(survival_pred.squeeze(), time)
            classification_loss = self.classification_criterion(classification_pred, grade)
            
            survival_losses.append(survival_loss.item())
            classification_losses.append(classification_loss.item())
            accuracies.append((classification_pred.argmax(dim=1) == grade).float().mean().item())
            
            survival_preds.extend(survival_pred.squeeze().cpu().numpy())
            survival_true.extend(time.cpu().numpy())
        
        c_index = concordance_index(survival_true, survival_preds)
        avg_survival_loss = sum(survival_losses) / len(survival_losses)
        avg_classification_loss = sum(classification_losses) / len(classification_losses)
        avg_accuracy = sum(accuracies) / len(accuracies)
        
        return {
            'val_survival_loss': avg_survival_loss,
            'val_classification_loss': avg_classification_loss,
            'val_accuracy': avg_accuracy,
            'val_c_index': c_index
        }
    
    def train(self, dataloaders: dict, epochs: int = 10, test_interval: int = 5):
        best_val_loss = float('inf')
        for epoch in range(epochs):
            train_metrics = []
            for batch in tqdm(dataloaders['train'], desc=f'Training Epoch {epoch+1}'):
                metrics = self.train_step(batch)
                train_metrics.append(metrics)
            
            avg_train_loss = sum(m['total_loss'] for m in train_metrics) / len(train_metrics)
            val_metrics = self.validate(dataloaders['validation'])
            
            self.scheduler.step(val_metrics['val_survival_loss'] + val_metrics['val_classification_loss'])
            
            print(f"Epoch {epoch+1}: Train Loss={avg_train_loss:.4f}, Val Loss={val_metrics['val_survival_loss'] + val_metrics['val_classification_loss']:.4f}, Accuracy={val_metrics['val_accuracy']:.4f}, C-Index={val_metrics['val_c_index']:.4f}")
            
            if val_metrics['val_survival_loss'] + val_metrics['val_classification_loss'] < best_val_loss:
                best_val_loss = val_metrics['val_survival_loss'] + val_metrics['val_classification_loss']
                torch.save(self.model.state_dict(), 'best_model.pt')

if __name__ == '__main__':
    from networks import RNASeqNet, load_data
    
    dataloaders = load_data('/Users/stanleychen/git/Melanoma/MultiCoFusion/data/processed_data.pkl')
    model = RNASeqNet(input_dim=10000)
    trainer = Trainer(model)
    trainer.train(dataloaders, epochs=50, test_interval=5)
