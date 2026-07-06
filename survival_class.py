import pandas as pd
import numpy as np

def add_survival_class(input_path, output_path):
    """
    Add survival_class column to patient data based on median survival time.
    
    Args:
        input_path (str): Path to input CSV file
        output_path (str): Path to save output CSV file
    """
    # Read the input CSV
    df = pd.read_csv(input_path)
    
    # Calculate median survival time
    median_survival = df['survival_time_days'].median()
    
    # Create new survival_class column (0 if <= median, 1 if > median)
    df['survival_class'] = (df['survival_time_days'] > median_survival).astype(int)
    
    # Save to new CSV file
    df.to_csv(output_path, index=False)
    
    # Print summary statistics
    print(f"Median survival time: {median_survival:.1f} days")
    print("\nClass distribution:")
    print(df['survival_class'].value_counts())
    
# Example usage
if __name__ == "__main__":
    input_path = "/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv"
    output_path = "/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv"
    
    add_survival_class(input_path, output_path)