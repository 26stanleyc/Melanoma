from torch.utils.data import DataLoader
import pickle
from data_loader import RNASeqDataset

def load_data(data_path):
    with open(data_path, 'rb') as f:
        data = pickle.load(f)
    
    datasets = {
        'train': DataLoader(RNASeqDataset(data['train']), batch_size=32, shuffle=True),
        'validation': DataLoader(RNASeqDataset(data['validation']), batch_size=32, shuffle=False),
        'test': DataLoader(RNASeqDataset(data['test']), batch_size=32, shuffle=False)
    }
    return datasets

data_path='/Users/stanleychen/git/Melanoma/MultiCoFusion/data/processed_data.pkl'
# Test DataLoader
dataloaders = load_data(data_path)

for batch in dataloaders['train']:
    print(batch['x_omic'].shape)
    break
