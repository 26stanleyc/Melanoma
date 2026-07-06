import torch
from torch.utils.data import Dataset, DataLoader
import pickle

from torch.utils.data import Dataset
import torch

from torch.utils.data import Dataset
import torch

class RNASeqDataset(Dataset):
    def __init__(self, data):
        """
        Custom Dataset for RNA-Seq Data.
        """
        self.x_omic = torch.tensor(data['x_omic_10000'], dtype=torch.float32)
        self.event = torch.tensor(data['e'], dtype=torch.float32)
        self.time = torch.tensor(data['t'], dtype=torch.float32)
        
        # Encode grades if they are strings
        if isinstance(data['g'][0], str):
            self.grade = torch.tensor(self.encode_grades(data['g']), dtype=torch.long)
        else:
            self.grade = torch.tensor(data['g'], dtype=torch.long)

    def __len__(self):
        return len(self.x_omic)

    def __getitem__(self, idx):
        if idx >= len(self):
            raise IndexError(f"Index {idx} out of range for dataset of length {len(self)}")
        return {
            'x_omic': self.x_omic[idx],
            'event': self.event[idx],
            'time': self.time[idx],
            'grade': self.grade[idx]
        }

    @staticmethod
    def encode_grades(grades):
        """
        Encode tumor grades as integers.
        """
        grade_mapping = {
            'Stage I': 0,
            'Stage II': 1,
            'Stage III': 2,
            'Stage IV': 3
        }
        return [grade_mapping.get(g, -1) for g in grades]



def load_data(data_path):
    """
    Load and preprocess RNA-Seq data from a pickle file.
    Args:
        data_path (str): Path to the pickle file containing data.
    Returns:
        dict: DataLoader objects for train, validation, and test datasets.
    """
    with open(data_path, 'rb') as f:
        data = pickle.load(f)

    datasets = {
        'train': RNASeqDataset(data['train']),
        'validation': RNASeqDataset(data['validation']),
        'test': RNASeqDataset(data['test'])
    }

    dataloaders = {
        phase: DataLoader(dataset, batch_size=32, shuffle=(phase == 'train'), drop_last=True)
        for phase, dataset in datasets.items()
    }

    return dataloaders

if __name__ == '__main__':
    data_path = '/Users/stanleychen/git/Melanoma/MultiCoFusion/data/processed_data.pkl'
    dataloaders = load_data(data_path)
    
    # Example usage
    for batch in dataloaders['train']:
        print(batch['x_omic'].shape, batch['event'].shape, batch['time'].shape, batch['grade'].shape)
        break
