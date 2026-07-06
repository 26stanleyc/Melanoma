import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.decomposition import PCA
import os

# Define your input and output paths here
INPUT_DATA_PATH = "/Users/stanleychen/git/Melanoma/apaper/miRNA.csv"
SURVIVAL_PATH = "/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv"
OUTPUT_DIR = "/Users/stanleychen/git/Melanoma/apaper"

def preprocess_survival_data(survival_path):
    """
    Preprocess survival data from the CSV file.
    """
    survival_df = pd.read_csv(survival_path)
    
    # Convert survival times to numeric
    survival_time = pd.to_numeric(survival_df['survival_time_days'])
    
    # Use event column directly
    survival_event = pd.to_numeric(survival_df['event'])
    
    # Set index to patient barcode
    survival_time.index = survival_df['patient_barcode']
    survival_event.index = survival_df['patient_barcode']
    
    print(f"Survival data loaded:")
    print(f"Number of patients: {len(survival_time)}")
    print(f"Number of events: {sum(survival_event)}")
    print(f"Mean survival time: {survival_time.mean():.1f} days")
    
    return survival_time, survival_event

def get_patient_id(barcode):
    """Extract patient ID from TCGA barcode (first 12 characters)"""
    return barcode[:12]

def calculate_cindex(risk_groups, survival_time, survival_event):
    """Calculate concordance index"""
    risk_scores = -1 * risk_groups
    c_index = concordance_index(survival_time, 
                              risk_scores,
                              survival_event)
    return c_index

def check_nans(data, stage_name):
    """Helper function to check for NaNs"""
    if isinstance(data, pd.DataFrame):
        nan_count = data.isna().sum().sum()
        if nan_count > 0:
            print(f"NaNs found at {stage_name}: {nan_count}")
            nan_cols = data.columns[data.isna().any()].tolist()
            print(f"Columns with NaNs: {nan_cols}")
    elif isinstance(data, np.ndarray):
        nan_count = np.isnan(data).sum()
        if nan_count > 0:
            print(f"NaNs found at {stage_name}: {nan_count}")
    elif isinstance(data, pd.Series):
        nan_count = data.isna().sum()
        if nan_count > 0:
            print(f"NaNs found at {stage_name}: {nan_count}")
    return nan_count > 0

def main():
    # Read data
    print("Reading data...")
    data = pd.read_csv(INPUT_DATA_PATH, index_col=0)
    print("\nInitial data structure:")
    print("Shape:", data.shape)
    print("Index (first 5):", data.index[:5].tolist())
    print("Columns (first 5):", data.columns[:5].tolist())
    
    # Always transpose to have samples as rows
    print("\nTransposing data to have samples as rows...")
    data = data.T
    print("New shape after transposition:", data.shape)
    print("New index (first 5):", data.index[:5].tolist())
    print("New columns (first 5):", data.columns[:5].tolist())
    
    # Read survival data
    survival_time, survival_event = preprocess_survival_data(SURVIVAL_PATH)
    
    # Extract patient IDs and match samples
    mirna_patient_ids = pd.Series(data.index).apply(get_patient_id)
    survival_patient_ids = pd.Series(survival_time.index).apply(get_patient_id)

    # Find common patients
    common_patients = set(mirna_patient_ids).intersection(set(survival_patient_ids))
    print(f"Number of common patients: {len(common_patients)}")

    # Create mapping from patient ID to full barcode for both datasets
    mirna_id_to_barcode = dict(zip(mirna_patient_ids, data.index))
    survival_id_to_barcode = dict(zip(survival_patient_ids, survival_time.index))

    # Get full barcodes for common patients
    mirna_barcodes = [mirna_id_to_barcode[pid] for pid in common_patients]
    survival_barcodes = [survival_id_to_barcode[pid] for pid in common_patients]

    # Filter data to keep only matched patients
    data = data.loc[mirna_barcodes]
    survival_time = survival_time[survival_barcodes]
    survival_event = survival_event[survival_barcodes]

    print(f"After matching - miRNA data shape: {data.shape}")
    print(f"After matching - number of survival times: {len(survival_time)}")
    print(f"After matching - number of survival events: {len(survival_event)}")

    # Scale data
    print("\nScaling data...")
    scaled_data = StandardScaler().fit_transform(data)
    print("Shape after scaling:", scaled_data.shape)
    print("NaNs after scaling:", np.isnan(scaled_data).sum())

    # Perform PCA
    print("\nPerforming PCA...")
    pca = PCA(n_components=210)
    reduced_data = pca.fit_transform(scaled_data)
    print("Shape after PCA:", reduced_data.shape)
    print("NaNs after PCA:", np.isnan(reduced_data).sum())

    # Convert to DataFrame
    reduced_df = pd.DataFrame(
        reduced_data,
        index=data.index,
        columns=[f'PC_{i+1}' for i in range(210)]
    )

    # Print variance explained
    print(f"\nTotal variance explained: {pca.explained_variance_ratio_.sum():.3f}")
    components_90 = np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.9) + 1
    print(f"Components needed for 90% variance: {components_90}")

    # Perform Cox-PH for each feature
    print("\nPerforming Cox-PH analysis...")
    significant_features = []
    p_values = []
    hazard_ratios = []

    for i, column in enumerate(reduced_df.columns):
        if i % 20 == 0:
            print(f"Processing feature {i+1}/{len(reduced_df.columns)}")
        
        try:
            cox_df = pd.DataFrame({
                'duration': survival_time,
                'event': survival_event,
                'feature': reduced_df[column]
            })
            
            cph = CoxPHFitter()
            cph.fit(cox_df, 'duration', 'event')
            
            # Access results directly from the model
            p_value = cph.print_summary(model_summary=True).loc['feature', 'p']
            hazard_ratio = np.exp(cph.params_['feature'])
            
            if p_value < 0.05:
                significant_features.append(column)
                p_values.append(p_value)
                hazard_ratios.append(hazard_ratio)
        except Exception as e:
            print(f"Error processing {column}: {str(e)}")
            continue

    print(f"\nFound {len(significant_features)} significant features")
    
    if len(significant_features) == 0:
        print("No significant features found. Try adjusting the p-value threshold.")
        return

    # Perform clustering
    print("\nPerforming clustering...")
    significant_data = reduced_df[significant_features]
    kmeans = KMeans(n_clusters=2, random_state=42)
    risk_groups = kmeans.fit_predict(significant_data)
    
    # Calculate C-index
    c_index = calculate_cindex(risk_groups, survival_time, survival_event)
    
    # Prepare results
    results_df = pd.DataFrame({
        'Sample': survival_barcodes,
        'Patient_ID': [get_patient_id(x) for x in survival_barcodes],
        'Risk_Group': ['High' if x == 1 else 'Low' for x in risk_groups],
        'Survival_Time': survival_time,
        'Event': survival_event
    })
    
    cox_results = pd.DataFrame({
        'Feature': significant_features,
        'P_Value': p_values,
        'Hazard_Ratio': hazard_ratios
    })

    # Save results
    print("\nSaving results...")
    results_df.to_csv(OUTPUT_DIR + '/risk_groups.csv')
    cox_results.to_csv(OUTPUT_DIR + '/cox_results.csv')
    
    # Save comprehensive metrics
    with open(OUTPUT_DIR + '/model_metrics.txt', 'w') as f:
        f.write(f"Concordance Index: {c_index:.3f}\n")
        f.write(f"Number of patients: {len(common_patients)}\n")
        f.write(f"Number of significant features: {len(significant_features)}\n")
        f.write(f"High risk patients: {sum(risk_groups == 1)}\n")
        f.write(f"Low risk patients: {sum(risk_groups == 0)}\n")
        f.write(f"\nPCA Metrics:\n")
        f.write(f"Total variance explained: {pca.explained_variance_ratio_.sum():.3f}\n")
        f.write(f"Components for 90% variance: {components_90}\n")
        
        # Add survival information
        f.write(f"\nSurvival Metrics:\n")
        f.write(f"Total events (deaths): {sum(survival_event)}\n")
        f.write(f"Mean survival time: {survival_time.mean():.1f} days\n")
    
    print("\nAnalysis Complete!")
    print(f"Found {len(significant_features)} significant features")
    print(f"High risk group: {sum(risk_groups == 1)} patients")
    print(f"Low risk group: {sum(risk_groups == 0)} patients")
    print(f"Concordance Index: {c_index:.3f}")

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    main()