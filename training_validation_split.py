import pandas as pd
import numpy as np

def create_train_val_split(input_file, output_file, train_size=260, random_seed=42):
    """
    Create training/validation split from patient data
    """
    # Set random seed for reproducibility
    np.random.seed(random_seed)
    
    # Read the data and handle quotes
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Extract patient IDs (skip empty lines)
    patient_ids = []
    for line in lines:
        parts = line.strip().split(',')
        if parts and parts[0].strip('"'):  # Only add non-empty IDs
            patient_ids.append(parts[0].strip('"'))
    
    # Random shuffle indices
    indices = np.arange(len(patient_ids))
    np.random.shuffle(indices)
    
    # Split into training and validation
    train_indices = indices[:train_size]
    val_indices = indices[train_size:]
    
    # Create dataframe with split assignments
    split_df = pd.DataFrame({
        'patient_id': patient_ids,
        'split': ['train' if i in train_indices else 'validation' for i in range(len(patient_ids))]
    })
    
    # Save splits to file
    split_df.to_csv(output_file, index=False)
    
    # Print summary
    print(f"Total patients: {len(patient_ids)}")
    print(f"Training patients: {len(train_indices)}")
    print(f"Validation patients: {len(val_indices)}")
    
    # Print first few assignments as verification
    print("\nFirst few training patients:")
    print(split_df[split_df['split'] == 'train']['patient_id'].head().values)
    print("\nFirst few validation patients:")
    print(split_df[split_df['split'] == 'validation']['patient_id'].head().values)

if __name__ == "__main__":
    create_train_val_split(
        input_file="/Users/stanleychen/git/Melanoma/preprocessing_data/RNA-seq-p/plier_encoded_data.csv",  # Update path
        output_file="training_validation_split.csv",
        train_size=260,
        random_seed=42
    )