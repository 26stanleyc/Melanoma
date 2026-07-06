import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt

def prepare_data(pathogenicity_file, labels_file):
    # Read the data
    path_scores = pd.read_csv(pathogenicity_file, header=None, 
                            names=['patient_id', 'pathogenicity_score'])
    labels = pd.read_csv(labels_file)
    
    # Merge the data
    data = pd.merge(path_scores, labels, 
                   left_on='patient_id', 
                   right_on='patient_barcode')
    
    # Create feature matrix X and target y
    X = data[['pathogenicity_score']]
    y = data['survival_time_days']
    events = (data['vital_status'] == 'Dead').astype(int)
    
    return X, y, events, data

def train_and_evaluate(X, y, events):
    # Split the data
    X_train, X_test, y_train, y_test, events_train, events_test = train_test_split(
        X, y, events, test_size=0.2, random_state=42
    )
    
    # Train linear regression model
    model = LinearRegression()
    model.fit(X_train, y_train)
    
    # Make predictions
    y_pred_train = model.predict(X_train)
    y_pred_test = model.predict(X_test)
    
    # Calculate metrics
    metrics = {
        'train': {
            'r2': r2_score(y_train, y_pred_train),
            'rmse': np.sqrt(mean_squared_error(y_train, y_pred_train)),
            'c_index': concordance_index(y_train, y_pred_train, events_train)
        },
        'test': {
            'r2': r2_score(y_test, y_pred_test),
            'rmse': np.sqrt(mean_squared_error(y_test, y_pred_test)),
            'c_index': concordance_index(y_test, y_pred_test, events_test)
        }
    }
    
    return model, metrics, (X_train, X_test, y_train, y_test, y_pred_train, y_pred_test)

def plot_predictions(results, save_path='prediction_plots.png'):
    X_train, X_test, y_train, y_test, y_pred_train, y_pred_test = results
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Training set
    ax1.scatter(y_train, y_pred_train, alpha=0.5)
    ax1.plot([y_train.min(), y_train.max()], [y_train.min(), y_train.max()], 'r--', lw=2)
    ax1.set_xlabel('Actual Survival Time (days)')
    ax1.set_ylabel('Predicted Survival Time (days)')
    ax1.set_title('Training Set Predictions')
    
    # Test set
    ax2.scatter(y_test, y_pred_test, alpha=0.5)
    ax2.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    ax2.set_xlabel('Actual Survival Time (days)')
    ax2.set_ylabel('Predicted Survival Time (days)')
    ax2.set_title('Test Set Predictions')
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def main():
    # Load and prepare data
    X, y, events, data = prepare_data('/Users/stanleychen/git/Melanoma/maf_sums.csv', '/Users/stanleychen/git/Melanoma/final_data/final_label_v2.csv')
    
    # Train model and get metrics
    model, metrics, results = train_and_evaluate(X, y, events)
    
    # Print metrics
    print("\nTraining Set Metrics:")
    print(f"R² Score: {metrics['train']['r2']:.3f}")
    print(f"RMSE: {metrics['train']['rmse']:.3f}")
    print(f"C-index: {metrics['train']['c_index']:.3f}")
    
    print("\nTest Set Metrics:")
    print(f"R² Score: {metrics['test']['r2']:.3f}")
    print(f"RMSE: {metrics['test']['rmse']:.3f}")
    print(f"C-index: {metrics['test']['c_index']:.3f}")
    
    # Plot predictions
    plot_predictions(results)
    
    return model, metrics, results

if __name__ == "__main__":
    model, metrics, results = main()