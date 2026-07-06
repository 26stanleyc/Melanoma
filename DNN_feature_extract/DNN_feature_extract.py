import tensorflow as tf
from tensorflow import keras
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import pickle
from tqdm import tqdm

class RNASeqFeatureExtractor:
    def __init__(self, input_dim):
        self.input_dim = input_dim
        
    def build_autoencoder(self):
        """
        Build autoencoder following the CNN architecture from the paper
        but adapted for RNA-seq data with correct dimension handling
        """
        # Input layer
        inputs = keras.Input(shape=(self.input_dim,))
        
        # Reshape for 1D convolution (adding channel dimension)
        x = keras.layers.Reshape((self.input_dim, 1))(inputs)
        
        # First CNN block
        x = keras.layers.Conv1D(filters=64, kernel_size=3, 
                              activation='relu', padding='same')(x)
        x = keras.layers.Conv1D(filters=64, kernel_size=3, 
                              activation='relu', padding='same')(x)
        x = keras.layers.MaxPooling1D(pool_size=2)(x)
        x = keras.layers.Dropout(0.4)(x)
        
        # Second CNN block
        x = keras.layers.Conv1D(filters=128, kernel_size=3, 
                              activation='relu', padding='same')(x)
        x = keras.layers.Conv1D(filters=128, kernel_size=3, 
                              activation='relu', padding='same')(x)
        x = keras.layers.MaxPooling1D(pool_size=2)(x)
        x = keras.layers.Dropout(0.4)(x)
        
        # Calculate dimensions after pooling
        # After two MaxPooling1D(2), dimension is reduced by factor of 4
        conv_output_dim = self.input_dim // 4
        
        # Flatten layer
        x = keras.layers.Flatten()(x)
        
        # Dense layers for encoding
        encoded = keras.layers.Dense(64, activation='relu', name='features')(x)
        
        # Decoder part (mirror of encoder)
        # Calculate the right dimension for reshaping
        decoder_dim = conv_output_dim * 128  # This should match the flattened dimension
        x = keras.layers.Dense(decoder_dim, activation='relu')(encoded)
        x = keras.layers.Reshape((conv_output_dim, 128))(x)
        
        # Decoder CNN blocks (transposed convolutions)
        x = keras.layers.Conv1DTranspose(filters=128, kernel_size=3, 
                                       strides=1, padding='same', 
                                       activation='relu')(x)
        x = keras.layers.Conv1DTranspose(filters=64, kernel_size=3, 
                                       strides=2, padding='same', 
                                       activation='relu')(x)
        
        # Final reconstruction
        decoded = keras.layers.Conv1DTranspose(filters=1, kernel_size=3, 
                                             strides=2, padding='same')(x)
        decoded = keras.layers.Flatten()(decoded)
        decoded = keras.layers.Dense(self.input_dim, activation='sigmoid')(decoded)
        
        # Create models
        self.autoencoder = keras.Model(inputs, decoded)
        self.encoder = keras.Model(inputs, encoded)
        
        # Print model summary for debugging
        print("\nEncoder Summary:")
        self.encoder.summary()
        print("\nAutoencoder Summary:")
        self.autoencoder.summary()
        
        return self.autoencoder, self.encoder

class ProgressCallback(keras.callbacks.Callback):
    def on_epoch_begin(self, epoch, logs=None):
        print(f'\nEpoch {epoch+1}/{self.params["epochs"]}')
        self.progbar = tqdm(total=self.params['steps'],
                          desc=f'Training',
                          unit='batch')

    def on_batch_end(self, batch, logs=None):
        self.progbar.update(1)
        self.progbar.set_postfix(loss=f'{logs["loss"]:.4f}')

    def on_epoch_end(self, epoch, logs=None):
        self.progbar.close()
        
def extract_features(expression_data, output_dir):
    """
    Extract features from RNA-seq data using the CNN autoencoder
    """
    print("Normalizing data...")
    scaler = StandardScaler()
    normalized_data = scaler.fit_transform(expression_data)
    
    print("Building model...")
    model = RNASeqFeatureExtractor(input_dim=expression_data.shape[1])
    autoencoder, encoder = model.build_autoencoder()
    
    # Compile
    autoencoder.compile(
        optimizer=keras.optimizers.Adam(),
        loss='mse'
    )
    
    # Train
    print("\nStarting training...")
    history = autoencoder.fit(
        normalized_data, 
        normalized_data,
        epochs=50,
        batch_size=32,
        validation_split=0.2,
        callbacks=[
            ProgressCallback(),
            keras.callbacks.EarlyStopping(
                monitor='val_loss',
                patience=5,
                restore_best_weights=True
            )
        ]
    )
    
    print("\nExtracting features...")
    features = encoder.predict(normalized_data)
    
    print("Saving results...")
    results = {
        'features': features,
        'feature_names': [f'feature_{i}' for i in range(features.shape[1])],
        'sample_ids': expression_data.index if isinstance(expression_data, pd.DataFrame) else None,
        'scaler': scaler,
        'model_architecture': encoder.get_config(),
        'training_history': history.history,
        'sample_info': {
            'total_samples': expression_data.shape[0],
            'feature_dimension': features.shape[1]
        }
    }
    
    # Save to files
    np.save(f'{output_dir}/features.npy', features)
    with open(f'{output_dir}/results.pkl', 'wb') as f:
        pickle.dump(results, f)
        
    print(f"Feature extraction complete! Extracted {features.shape[1]} features from {features.shape[0]} samples.")
    return results



if __name__ == "__main__":
    print("=== Starting RNA-Seq Feature Extraction Pipeline ===\n")
    
    # Define paths
    tcga_path = '/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/filtered_rna_seq-R.csv'
    output_dir = '/Users/stanleychen/git/Melanoma/DNN_feature_extract'
    
    # Load data
    print("Loading data...")
    expression_data = pd.read_csv(tcga_path, sep=',', index_col=0)
    expression_data = expression_data.T  # Transpose after loading
    
    # Extract features
    print("\n=== Starting Feature Extraction ===\n")
    results = extract_features(expression_data, output_dir)
    
    print("\n=== Pipeline Complete ===")
    print(f"Results saved to: {output_dir}")
    print(f"Features shape: {results['features'].shape}")
# Access features
#features = results['features']
