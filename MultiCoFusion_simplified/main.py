import argparse
import logging
from pathlib import Path
import yaml
import torch
from torch.utils.data import DataLoader
import pandas as pd
import numpy as np
from datetime import datetime

from model import FlexibleOmicsModel
from dataset import OmicsDataset
from training import Trainer
from validation import ModelValidator


def setup_logging(output_dir: Path):
    """Setup logging configuration."""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = output_dir / f'run_{timestamp}.log'
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Add timestamp to output directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    config['output_dir'] = str(Path(config['output_dir']) / timestamp)
    
    return config

def load_data(config: dict) -> tuple:
    """Load and prepare data according to configuration."""
    logging.info("Loading clinical data...")
    # Load clinical data
    clinical_df = pd.read_csv(config['data']['clinical_path'])
    
    # Replace hyphens with dots in patient barcodes to match RNA-seq data format
    clinical_df['patient_barcode'] = clinical_df['patient_barcode'].str.replace('-', '.')
    clinical_df.set_index('patient_barcode', inplace=True)
    
    # Filter clinical data
    clinical_df = clinical_df[
        (clinical_df['split'].notna()) &  # Has split information
        (clinical_df['tumor_stage'] != '[Not Available]') &
        (clinical_df['survival_time_days'] != '[Not Available]') &
        (~clinical_df['tumor_stage'].str.contains('Stage 0')) &  # Exclude Stage 0
        (~clinical_df['tumor_stage'].str.contains('Stage0'))     # Handle potential variations
    ]
    
    logging.info(f"Clinical data loaded with {len(clinical_df)} samples after filtering")
    
    # Load omic data
    logging.info("Loading omic data...")
    omic_data = {}
    for omic_type, path in config['data']['omic_paths'].items():
        if path:  # Only load if path is provided
            logging.info(f"Loading {omic_type} data from {path}")
            df = pd.read_csv(path, index_col=0)
            
            # Find common indices
            common_indices = clinical_df.index.intersection(df.index)
            if len(common_indices) == 0:
                raise ValueError(f"No common samples found between clinical data and {omic_type} data")
            
            logging.info(f"Found {len(common_indices)} common samples for {omic_type}")
            
            # Filter both dataframes to common indices
            clinical_df = clinical_df.loc[common_indices]
            df = df.loc[common_indices]
            
            omic_data[omic_type] = df
    
    if not omic_data:
        raise ValueError("No omic data was loaded. Check your file paths.")
    
    logging.info(f"Final dataset size: {len(clinical_df)} samples")
    return clinical_df, omic_data

def save_experiment_config(config: dict, output_dir: Path):
    """Save experiment configuration."""
    with open(output_dir / 'experiment_config.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

def main():
    parser = argparse.ArgumentParser(description='Train and validate multi-omics model')
    parser.add_argument('--config', type=str, required=True, help='Path to config file')
    parser.add_argument('--mode', choices=['train', 'validate', 'predict'], required=True)
    parser.add_argument('--model_path', type=str, help='Path to model for validation/prediction')
    args = parser.parse_args()

    # Load configuration
    config = load_config(args.config)
    
    # Setup output directory
    output_dir = Path(config['output_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save experiment configuration
    save_experiment_config(config, output_dir)
    
    # Setup logging
    setup_logging(output_dir)
    logger = logging.getLogger(__name__)
    
    # Set random seeds for reproducibility
    torch.manual_seed(config['seed'])
    np.random.seed(config['seed'])
    
    # Load data
    logger.info("Loading data...")
    clinical_df, omic_data = load_data(config)
    
    # Create dataset
    dataset = OmicsDataset(
        clinical_data=clinical_df,
        omic_data=omic_data,
        mode=args.mode
    )
    
    # Create model
    input_dims = {
        omic_type: data.shape[1] 
        for omic_type, data in omic_data.items()
    }
    
    model = FlexibleOmicsModel(
        modality_dims=config['model']['modality_dims'],
        fusion_dim=config['model']['fusion_dim'],
        num_classes=config['model']['num_classes'],
        dropout=config['model']['dropout']
    )
    
    if args.mode == 'train':
        # Training mode
        trainer = Trainer(
            model=model,
            lr=config['training']['learning_rate'],
            weight_decay=config['training']['weight_decay'],
            alpha=config['training']['alpha']
        )
        
        # Create data loader
        train_loader = DataLoader(
            dataset,
            batch_size=config['training']['batch_size'],
            shuffle=True,
            num_workers=config['training']['num_workers'],
            pin_memory=True
        )
        
        # Train with cross-validation
        logger.info("Starting training with cross-validation...")
        cv_results = trainer.train_with_cross_validation(
            dataset=dataset,
            save_dir=output_dir / 'models',
            n_splits=config['training']['n_splits'],
            num_epochs=config['training']['num_epochs'],
            batch_size=config['training']['batch_size'],
            early_stopping_patience=config['training']['early_stopping_patience']
        )
        
        logger.info("Training completed successfully.")
        
    elif args.mode == 'validate':
        if not args.model_path:
            raise ValueError("Model path must be provided for validation mode")
        
        validator = ModelValidator(model=model)
        
        # Create data loader
        val_loader = DataLoader(
            dataset,
            batch_size=config['validation']['batch_size'],
            shuffle=False,
            num_workers=config['validation']['num_workers'],
            pin_memory=True
        )
        
        # Validate model
        logger.info("Starting validation...")
        metrics = validator.validate_model(
            val_loader=val_loader,
            model_dir=args.model_path,
            output_dir=output_dir / 'validation'
        )
        
        logger.info("Validation completed successfully.")
        
    else:  # predict mode
        if not args.model_path:
            raise ValueError("Model path must be provided for prediction mode")
        
        validator = ModelValidator(model=model)
        
        # Create data loader
        test_loader = DataLoader(
            dataset,
            batch_size=config['prediction']['batch_size'],
            shuffle=False,
            num_workers=config['prediction']['num_workers'],
            pin_memory=True
        )
        
        # Make predictions
        logger.info("Starting prediction...")
        stage_preds, survival_preds = validator.predict(
            test_loader=test_loader,
            model_path=args.model_path
        )
        
        # Save predictions
        predictions_df = pd.DataFrame({
            'patient_id': clinical_df.index,
            'predicted_stage': stage_preds,
            'predicted_survival': survival_preds.flatten()
        })
        predictions_df.to_csv(output_dir / 'predictions.csv', index=False)
        
        logger.info("Prediction completed successfully.")

if __name__ == '__main__':
    main()