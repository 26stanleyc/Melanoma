import torch
import logging
from pathlib import Path
from networks import RNASeqNet, load_data
from training import Trainer
from validate import validate_model

# Setup Logging
def setup_logging(log_dir: str):
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=Path(log_dir) / 'training.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    logging.getLogger('').addHandler(console)

# Main Script
def main():
    # Paths
    data_path = '/Users/stanleychen/git/Melanoma/MultiCoFusion/data/processed_data.pkl'
    model_save_dir = './models'
    log_dir = './logs'
    
    Path(model_save_dir).mkdir(parents=True, exist_ok=True)
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    
    # Setup Logging
    setup_logging(log_dir)
    logging.info('Starting RNA-Seq Training Pipeline')
    
    # Load Data
    logging.info('Loading data...')
    dataloaders = load_data(data_path)
    
    # Initialize Model
    logging.info('Initializing model...')
    model = RNASeqNet(input_dim=10000)
    
    # Initialize Trainer
    trainer = Trainer(model)
    
    # Train the Model
    logging.info('Starting training...')
    trainer.train(dataloaders, epochs=50, test_interval=5)
    
    # Save Final Model
    final_model_path = Path(model_save_dir) / 'final_model.pt'
    torch.save(model.state_dict(), final_model_path)
    logging.info(f'Training completed. Final model saved at {final_model_path}')
    
    # Validate the Model
    logging.info('Starting validation...')
    validation_results = validate_model(model, dataloaders['validation'])
    
    # Save Validation Metrics
    metrics_path = Path(log_dir) / 'validation_metrics.txt'
    with open(metrics_path, 'w') as f:
        for key, value in validation_results.items():
            f.write(f'{key}: {value:.4f}\n')
    logging.info(f'Validation completed. Metrics saved at {metrics_path}')
    
    logging.info('Pipeline completed successfully.')

if __name__ == '__main__':
    main()
