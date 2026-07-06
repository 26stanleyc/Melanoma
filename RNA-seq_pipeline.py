import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from torchvision import transforms
from PIL import Image
from transformers import ResNetModel, ResNetConfig
import pandas as pd
import os

# Define directories and parameters
image_data_dir = "/Users/stanleychen/git/Melanoma/final_data/image_data"
plier_data_file = "/Users/stanleychen/git/Melanoma/preprocessing_data/RNA-seq-p/plier_encoded_data.csv"  # Path to the PLIER data file
label_file = "/Users/stanleychen/git/Melanoma/final_data/final_label_v2.csv"  # Path to the label file
batch_size = 32
learning_rate = 1e-4
epochs = 10

# Load and preprocess label data
def load_label_data(file_path):
    label_data = pd.read_csv(file_path)
    label_data = label_data[label_data["split"] == "train"]  # Filter for training data only
    label_data["patient_barcode"] = label_data["patient_barcode"].str.replace("-", ".")  # Normalize patient IDs
    return label_data

label_data = load_label_data(label_file)

# Custom Dataset to handle image-label mapping
class CustomImageDataset(Dataset):
    def __init__(self, image_dir, label_data, transform=None):
        self.image_dir = image_dir
        self.label_data = label_data
        self.transform = transform

        # Normalize patient IDs to ensure consistency
        self.image_paths = []
        self.labels = []
        unmatched_ids = []
        for _, row in label_data.iterrows():
            patient_id = row["patient_barcode"].replace(".", "-")  # Convert to match image filenames
            tumor_stage = row["tumor_stage"]
            survival_time = row["survival_time_days"]

            # Check for matching image files
            matched = False
            for file_name in os.listdir(image_dir):
                if patient_id in file_name:
                    self.image_paths.append(os.path.join(image_dir, file_name))
                    self.labels.append((tumor_stage, survival_time))
                    matched = True
                    break
            if not matched:
                unmatched_ids.append(patient_id)

        if unmatched_ids:
            print(f"Warning: No matching images found for patient IDs: {unmatched_ids}")

    def __len__(self):
        return len(self.image_paths)

    def __getitem__(self, idx):
        image_path = self.image_paths[idx]
        image = Image.open(image_path).convert("RGB")
        label = self.labels[idx]

        if self.transform:
            image = self.transform(image)

        # Split the labels into tumor stage (categorical) and survival time (regression)
        tumor_stage, survival_time = label
        tumor_stage = torch.tensor(int(tumor_stage), dtype=torch.long)
        survival_time = torch.tensor(float(survival_time), dtype=torch.float32)

        return image, tumor_stage, survival_time



# Define image transformations
transform = transforms.Compose([
    transforms.Resize((224, 224)),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
])

# Create the custom dataset and dataloader
custom_dataset = CustomImageDataset(image_data_dir, label_data, transform=transform)
dataloader = DataLoader(custom_dataset, batch_size=batch_size, shuffle=True)

# Load and preprocess PLIER data
def load_plier_data(file_path, patient_ids):
    plier_data = pd.read_csv(file_path, index_col=0)
    plier_data.index = plier_data.index.str.replace(".", "-")  # Normalize PLIER IDs to match label file format
    plier_data = plier_data.loc[plier_data.index.intersection(patient_ids)]
    plier_data = torch.tensor(plier_data.values, dtype=torch.float32)
    return plier_data

plier_data_tensor = load_plier_data(plier_data_file, label_data["patient_barcode"].values)

# Load pre-trained ResNet-152 model from Hugging Face
resnet_model = ResNetModel.from_pretrained("microsoft/resnet-152")

# Freeze ResNet feature extraction layers
for param in resnet_model.parameters():
    param.requires_grad = False

# Define a fusion layer (FCNN) for embedding fusion
class FusionModel(nn.Module):
    def __init__(self, resnet_embedding_size, plier_embedding_size, fused_embedding_size, num_tumor_stages):
        super(FusionModel, self).__init__()
        self.resnet_model = resnet_model
        self.fc_fusion = nn.Sequential(
            nn.Linear(resnet_embedding_size + plier_embedding_size, 512),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(512, fused_embedding_size),
            nn.ReLU()
        )
        self.tumor_stage_head = nn.Linear(fused_embedding_size, num_tumor_stages)  # For tumor stage classification
        self.survival_time_head = nn.Linear(fused_embedding_size, 1)  # For survival time regression

    def forward(self, images, plier_embeddings):
        # Extract image embeddings from ResNet
        image_embeddings = self.resnet_model(images).pooler_output
        # Concatenate image embeddings and PLIER embeddings
        fused_input = torch.cat((image_embeddings, plier_embeddings), dim=1)
        # Pass through the fusion layer
        fused_output = self.fc_fusion(fused_input)
        # Predict tumor stage and survival time
        tumor_stage_output = self.tumor_stage_head(fused_output)
        survival_time_output = self.survival_time_head(fused_output).squeeze()
        return tumor_stage_output, survival_time_output

# Define model parameters
resnet_embedding_size = resnet_model.config.hidden_sizes[-1]
  # ResNet output embedding size from Hugging Face config
plier_embedding_size = plier_data_tensor.size(1)  # Dimensionality of PLIER embeddings
print(f"ResNet Embedding Size: {resnet_embedding_size}, Type: {type(resnet_embedding_size)}")
print(f"PLIER Embedding Size: {plier_embedding_size}, Type: {type(plier_embedding_size)}")

fused_embedding_size = 256
num_tumor_stages = len(set(label_data["tumor_stage"].astype("category").cat.codes))  # Number of unique tumor stages

# Instantiate the fusion model
model = FusionModel(resnet_embedding_size, plier_embedding_size, fused_embedding_size, num_tumor_stages)

# Define loss functions and optimizer
classification_loss_fn = nn.CrossEntropyLoss()
regression_loss_fn = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Fine-tune the model
for epoch in range(epochs):
    model.train()
    for batch_idx, (images, tumor_stage_labels, survival_time_labels) in enumerate(dataloader):
        # Prepare PLIER data and labels for the batch
        batch_plier_data = plier_data_tensor[batch_idx * batch_size:(batch_idx + 1) * batch_size].cuda()
        tumor_stage_labels = tumor_stage_labels.cuda()
        survival_time_labels = survival_time_labels.cuda()

        # Move data to GPU if available
        images = images.cuda()
        
        # Forward pass
        tumor_stage_output, survival_time_output = model(images, batch_plier_data)

        # Compute losses
        classification_loss = classification_loss_fn(tumor_stage_output, tumor_stage_labels)
        regression_loss = regression_loss_fn(survival_time_output, survival_time_labels)
        loss = classification_loss + regression_loss

        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    print(f"Epoch [{epoch + 1}/{epochs}], Loss: {loss.item():.4f}")

# Save the fine-tuned model
torch.save(model.state_dict(), "fine_tuned_resnet_plier_fusion_model.pth")

print("Model fine-tuning and fusion layer training completed.")
