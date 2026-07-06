import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from lifelines.utils import concordance_index
import pickle
from torch.utils.data import Dataset, DataLoader

class GSNN(nn.Module):
    def __init__(self, input_dim=10000, omic_dim=32, dropout_rate=0.25, act_1=None, act_2=None, label_dim_1=1, label_dim_2=2, init_max=True):
        super(GSNN, self).__init__()
        hidden = [64, 48, 32, 32]
        self.act_1 = act_1 if act_1 else nn.ELU()
        self.act_2 = act_2 if act_2 else nn.ELU()

        self.encoder0 = nn.Linear(input_dim, input_dim)
        self.encoder0_elu = nn.ELU()
        self.encoder0_alphadropout = nn.AlphaDropout(p=dropout_rate, inplace=False)

        encoder1 = nn.Sequential(
            nn.Linear(input_dim, hidden[0]),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder2 = nn.Sequential(
            nn.Linear(hidden[0], hidden[1]),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder3 = nn.Sequential(
            nn.Linear(hidden[1], hidden[2]),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder4 = nn.Sequential(
            nn.Linear(hidden[2], omic_dim),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder1_2 = nn.Sequential(
            nn.Linear(input_dim, hidden[0]),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder2_2 = nn.Sequential(
            nn.Linear(hidden[0], hidden[1]),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder3_2 = nn.Sequential(
            nn.Linear(hidden[1], hidden[2]),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        encoder4_2 = nn.Sequential(
            nn.Linear(hidden[2], omic_dim),
            nn.ELU(),
            nn.AlphaDropout(p=dropout_rate, inplace=False))

        self.encoder_1 = nn.Sequential(encoder1, encoder2, encoder3, encoder4)
        self.encoder_2 = nn.Sequential(encoder1_2, encoder2_2, encoder3_2, encoder4_2)
        self.classifier_1 = nn.Sequential(nn.Linear(omic_dim, label_dim_1))
        self.classifier_2 = nn.Sequential(nn.Linear(omic_dim, label_dim_2))

        if init_max:
            self.init_max_weights()

        self.output_range = nn.Parameter(torch.FloatTensor([6]), requires_grad=False)
        self.output_shift = nn.Parameter(torch.FloatTensor([-3]), requires_grad=False)

    def init_max_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(self, x):
        x = self.encoder0(x)
        x = self.encoder0_elu(x)
        x = self.encoder0_alphadropout(x)

        x1 = self.encoder_1(x)
        x2 = self.encoder_2(x)
        x = x1 + x2

        survival_pred = self.classifier_1(x)
        grade_pred = self.classifier_2(x)
        
        return {
            'survival': survival_pred,
            'grade': grade_pred
        }

class SGCN(nn.Module):
    def __init__(self, hidden_dim, output_dim, dropout=0.3):
        super(SGCN, self).__init__()
        self.conv1 = GCNConv(hidden_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, output_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, edge_index=None):
        if edge_index is None:
            edge_index = torch.tensor([[i, i+1] for i in range(x.size(0)-1)], dtype=torch.long).t().contiguous()
        
        x = F.relu(self.conv1(x, edge_index))
        x = self.dropout(x)
        x = F.relu(self.conv2(x, edge_index))
        return x

class RNASeqNet(nn.Module):
    def __init__(self, input_dim, hidden_dim=128, output_dim=64, dropout=0.3):
        super(RNASeqNet, self).__init__()
        self.gsnn = GSNN(input_dim=input_dim, omic_dim=hidden_dim, dropout_rate=dropout)
        self.sgcn = SGCN(hidden_dim=hidden_dim, output_dim=output_dim, dropout=dropout)
        self.survival_head = nn.Linear(output_dim, 1)
        self.grade_head = nn.Linear(output_dim, 3)

    def forward(self, x_omic):
        x = self.gsnn(x_omic)['survival']
        x = self.sgcn(x)
        survival_pred = self.survival_head(x)
        grade_pred = self.grade_head(x)
        
        return {
            'survival': survival_pred,
            'grade': grade_pred
        }

def calculate_metrics(y_true, y_pred):
    acc = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average='macro')
    recall = recall_score(y_true, y_pred, average='macro')
    f1 = f1_score(y_true, y_pred, average='macro')
    return acc, precision, recall, f1

def load_data(data_path):
    with open(data_path, 'rb') as f:
        data = pickle.load(f)
    datasets = {
        'train': DataLoader(data['train'], batch_size=32, shuffle=True),
        'validation': DataLoader(data['validation'], batch_size=32, shuffle=False),
        'test': DataLoader(data['test'], batch_size=32, shuffle=False)
    }
    return datasets
