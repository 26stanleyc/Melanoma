import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.decomposition import PCA

# Define your input and output paths here
INPUT_DATA_PATH = "/Users/stanleychen/git/Melanoma/apaper/miRNA.csv"
SURVIVAL_TIME_PATH = "/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv"
SURVIVAL_EVENT_PATH = "/Users/stanleychen/git/Melanoma/afpipeline/patient_survival_obs.csv"
OUTPUT_DIR = "/Users/stanleychen/git/Melanoma/apaper"

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

def main():
    # Read data
    print("Reading data...")
    data = pd.read_csv(INPUT_DATA_PATH, index_col=0)
    survival_time = pd.read_csv(SURVIVAL_TIME_PATH, index_col=0).iloc[:, 0]
    survival_event = pd.read_csv(SURVIVAL_EVENT_PATH, index_col=0).iloc[:, 0]

    # Check initial data for NaNs
    check_nans(data, "Initial miRNA data")
    check_nans(survival_time, "Initial survival time")
    check_nans(survival_event, "Initial survival event")

    print(f"Initial miRNA data shape: {data.shape}")
    print(f"Initial survival time shape: {survival_time.shape}")
    print(f"Initial survival event shape: {survival_event.shape}")

    # Extract patient IDs and match samples
    mirna_patient_ids = pd.Series(data.columns).apply(get_patient_id)
    survival_patient_ids = pd.Series(survival_time.index).apply(get_patient_id)

    common_patients = set(mirna_patient_ids).intersection(set(survival_patient_ids))
    print(f"Number of common patients: {len(common_patients)}")

    mirna_id_to_barcode = dict(zip(mirna_patient_ids, data.columns))
    survival_id_to_barcode = dict(zip(survival_patient_ids, survival_time.index))

    mirna_barcodes = [mirna_id_to_barcode[pid] for pid in common_patients]
    survival_barcodes = [survival_id_to_barcode[pid] for pid in common_patients]

    # Filter data to keep only matched patients
    data = data[mirna_barcodes]
    survival_time = survival_time[survival_barcodes]
    survival_event = survival_event[survival_barcodes]

    # Check data after matching
    check_nans(data, "After matching - miRNA data")
    check_nans(survival_time, "After matching - survival time")
    check_nans(survival_event, "After matching - survival event")

    print(f"After matching - miRNA data shape: {data.shape}")
    print(f"After matching - number of survival times: {len(survival_time)}")
    print(f"After matching - number of survival events: {len(survival_event)}")

    # Verify matching
    print("\nVerifying patient matching...")
    matched_mirna_ids = pd.Series(data.columns).apply(get_patient_id)
    matched_survival_ids = pd.Series(survival_time.index).apply(get_patient_id)
    all_matched = all(matched_mirna_ids == matched_survival_ids)
    print(f"All patients properly matched: {all_matched}")

    # Scale and perform PCA
    print("\nPerforming scaling and PCA...")
    scaled_data = StandardScaler().fit_transform(data)
    check_nans(scaled_data, "After scaling")

    pca = PCA(n_components=210)
    reduced_data = pca.fit_transform(scaled_data)
    check_nans(reduced_data, "After PCA")

    reduced_df = pd.DataFrame(
        reduced_data, 
        columns=[f'PC_{i+1}' for i in range(210)], 
        index=data.index
    )
    check_nans(reduced_df, "After converting PCA to DataFrame")

    # Print PCA stats
    print("\nPCA statistics:")
    print(f"Explained variance ratio sum: {np.sum(pca.explained_variance_ratio_):.3f}")
    print(f"Number of components explaining 90% variance: {np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.9) + 1}")

    # Perform Cox-PH for each feature
    print("\nPerforming Cox-PH analysis...")
    significant_features = []
    p_values = []
    hazard_ratios = []

    for i, column in enumerate(reduced_df.columns):
        if i % 20 == 0:  # Print progress
            print(f"Processing feature {i+1}/{len(reduced_df.columns)}")
            
        try:
            cox_df = pd.DataFrame({
                'duration': survival_time,
                'event': survival_event,
                'feature': reduced_df[column]
            })
            
            if cox_df.isna().any().any():
                print(f"NaNs found in Cox data for {column}")
                continue

            cph = CoxPHFitter()
            cph.fit(cox_df, 'duration', 'event')
            summary = cph.print_summary()
            
            p_value = summary.iloc[0]['p']
            hazard_ratio = np.exp(cph.params_.iloc[0])
            
            if p_value < 0.05:
                significant_features.append(column)
                p_values.append(p_value)
                hazard_ratios.append(hazard_ratio)
        except Exception as e:
            print(f"Error processing {column}: {str(e)}")
            continue

    print(f"\nFound {len(significant_features)} significant features")
    
    # Check if we have any significant features
    if len(significant_features) == 0:
        print("No significant features found. Try adjusting the p-value threshold.")
        return
        
    # Select significant features and perform clustering
    significant_data = reduced_df[significant_features]
    print("\nPerforming clustering...")
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
    
    # Save metrics
    with open(OUTPUT_DIR + '/model_metrics.txt', 'w') as f:
        f.write(f"Concordance Index: {c_index:.3f}\n")
        f.write(f"Number of patients: {len(common_patients)}\n")
        f.write(f"Number of significant features: {len(significant_features)}\n")
        f.write(f"High risk patients: {sum(risk_groups == 1)}\n")
        f.write(f"Low risk patients: {sum(risk_groups == 0)}\n")
        f.write(f"\nPCA Metrics:\n")
        f.write(f"Total variance explained: {np.sum(pca.explained_variance_ratio_):.3f}\n")
        f.write(f"Components for 90% variance: {np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.9) + 1}\n")
    
    # Print summary
    print("\nAnalysis Complete!")
    print(f"Found {len(significant_features)} significant features")
    print(f"High risk group: {sum(risk_groups == 1)} patients")
    print(f"Low risk group: {sum(risk_groups == 0)} patients")
    print(f"Concordance Index: {c_index:.3f}")

if __name__ == "__main__":
    # Create output directory if it doesn't exist
    import os
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Run analysis
    main()