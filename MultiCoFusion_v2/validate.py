import torch
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from lifelines.utils import concordance_index
from networks import RNASeqNet, load_data

@torch.no_grad()
def validate_model(model: RNASeqNet, dataloader: DataLoader, device='cuda' if torch.cuda.is_available() else 'cpu'):
    model.to(device)
    model.eval()
    
    survival_losses = []
    classification_losses = []
    accuracies = []
    precisions = []
    recalls = []
    f1_scores = []
    survival_preds = []
    survival_true = []
    grade_true = []
    grade_preds = []
    
    survival_criterion = torch.nn.MSELoss(reduction='mean')
    classification_criterion = torch.nn.CrossEntropyLoss()
    
    for batch in dataloader:
        x_omic = batch['x_omic'].to(device)
        time = batch['time'].to(device)
        event = batch['event'].to(device)
        grade = batch['grade'].to(device)
        
        outputs = model(x_omic)
        survival_loss = survival_criterion(outputs['survival'].squeeze(), time)
        classification_loss = classification_criterion(outputs['grade'], grade)
        
        survival_losses.append(survival_loss.item())
        classification_losses.append(classification_loss.item())
        
        preds = outputs['grade'].argmax(dim=1)
        accuracies.append((preds == grade).float().mean().item())
        precisions.append(precision_score(grade.cpu(), preds.cpu(), average='macro', zero_division=0))
        recalls.append(recall_score(grade.cpu(), preds.cpu(), average='macro', zero_division=0))
        f1_scores.append(f1_score(grade.cpu(), preds.cpu(), average='macro', zero_division=0))
        
        survival_preds.extend(outputs['survival'].squeeze().cpu().numpy())
        survival_true.extend(time.cpu().numpy())
        grade_true.extend(grade.cpu().numpy())
        grade_preds.extend(preds.cpu().numpy())
    
    c_index = concordance_index(survival_true, survival_preds)
    avg_survival_loss = sum(survival_losses) / len(survival_losses)
    avg_classification_loss = sum(classification_losses) / len(classification_losses)
    avg_accuracy = sum(accuracies) / len(accuracies)
    avg_precision = sum(precisions) / len(precisions)
    avg_recall = sum(recalls) / len(recalls)
    avg_f1_score = sum(f1_scores) / len(f1_scores)
    
    results = {
        'validation_survival_loss': avg_survival_loss,
        'validation_classification_loss': avg_classification_loss,
        'validation_accuracy': avg_accuracy,
        'validation_precision': avg_precision,
        'validation_recall': avg_recall,
        'validation_f1_score': avg_f1_score,
        'validation_c_index': c_index
    }
    
    print("Validation Results:")
    for key, value in results.items():
        print(f"{key}: {value:.4f}")
    
    return results

if __name__ == '__main__':
    dataloaders = load_data('/Users/stanleychen/git/Melanoma/MultiCoFusion/data/processed_data.pkl')
    model = RNASeqNet(input_dim=10000)
    model.load_state_dict(torch.load('best_model.pt'))
    validate_model(model, dataloaders['validation'])
