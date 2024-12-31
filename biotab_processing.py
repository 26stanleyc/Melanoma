import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

def process_melanoma_data(file_path):
    # Load the data
    data = pd.read_csv(file_path, sep='\t', skiprows=3)

    # Define feature categories
    patient_history = [
        'birth_days_to', 'gender', 'height_cm_at_diagnosis', 'weight_kg_at_diagnosis',
        'race', 'ethnicity', 'history_other_malignancy', 'primary_melanoma_known_dx',
        'primary_multiple_at_dx', 'primary_at_dx_count', 'breslow_thickness_at_diagnosis',
        'clark_level_at_diagnosis', 'primary_melanoma_tumor_ulceration',
        'primary_melanoma_mitotic_rate', 'age_at_diagnosis', 'ldh_level',
        'ajcc_pathologic_tumor_stage', 'melanoma_primary_count', 'clinical_M',
        'clinical_N', 'clinical_T', 'clinical_stage'
    ]

    patient_treatment = [
        'history_neoadjuvant_treatment', 'radiation_therapy_to_primary',
        'prior_radiation_therapy', 'history_neoadjuvant_tx_type',
        'ifn_tx_90_days_prior_to_resection', 'radiation_treatment_adjuvant',
        'pharmaceutical_tx_adjuvant'
    ]

    labels = [
        'tumor_status', 'vital_status', 'new_tumor_event_dx_indicator',
        'days_to_patient_progression_free', 'days_to_tumor_progression'
    ]

    # Combine features
    features = patient_history + patient_treatment

    # Identify numeric and categorical columns
    numeric_features = data[features].select_dtypes(include=['int64', 'float64']).columns
    categorical_features = data[features].select_dtypes(include=['object']).columns

    # Create preprocessing pipelines
    numeric_transformer = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='median')),
        ('scaler', StandardScaler())
    ])

    categorical_transformer = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
        ('onehot', OneHotEncoder(handle_unknown='ignore'))
    ])

    # Combine preprocessing steps
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', numeric_transformer, numeric_features),
            ('cat', categorical_transformer, categorical_features)
        ])

    # Fit and transform the features
    X = preprocessor.fit_transform(data[features])

    # Process labels
    y = data[labels]

    # Handle categorical labels
    y['tumor_status'] = y['tumor_status'].map({'TUMOR FREE': 0, 'WITH TUMOR': 1})
    y['vital_status'] = y['vital_status'].map({'Alive': 0, 'Dead': 1})
    y['new_tumor_event_dx_indicator'] = y['new_tumor_event_dx_indicator'].map({'NO': 0, 'YES': 1})

    # Convert to numpy arrays
    y = y.to_numpy()

    return X, y, preprocessor

# Usage example:
X, y, preprocessor = process_melanoma_data('/Users/stanleychen/git/Melanoma/test_data/gdc_download_20240809_223437.267906/58cbbc07-5ec4-47c7-9295-11ccbf7693f4/nationwidechildrens.org_clinical_patient_skcm.txt')
print(y)