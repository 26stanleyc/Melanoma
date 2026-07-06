import torch
from torch.utils.data import Dataset
import pandas as pd
import numpy as np
from typing import Dict, Optional, List
from sklearn.preprocessing import StandardScaler, RobustScaler
import logging

class OmicsDataset(Dataset):
    def __init__(
        self,
        clinical_data: pd.DataFrame,
        omic_data: Dict[str, pd.DataFrame],
        mode: str = 'train',
        scalers: Optional[Dict[str, StandardScaler]] = None
    ):
        self.clinical_data = clinical_data
        self.mode = mode
        self.logger = logging.getLogger(__name__)
        
        # Normalize survival times
        survival_scaler = StandardScaler()
        if mode == 'train':
            self.survival_times = survival_scaler.fit_transform(
                clinical_data['survival_time_days'].values.reshape(-1, 1)
            ).flatten()
            self.survival_scaler = survival_scaler
        else:
            self.survival_times = survival_scaler.transform(
                clinical_data['survival_time_days'].values.reshape(-1, 1)
            ).flatten()
        
        # Process stages
        self.stage_labels = self._process_stages(clinical_data['tumor_stage'])
        
        # Process and scale omic data
        self.omic_data = {}
        self.scalers = {}
        
        for modality, data in omic_data.items():
            if scalers is None:
                scaler = RobustScaler()
                scaled_data = scaler.fit_transform(data)
                self.scalers[modality] = scaler
            else:
                scaler = scalers[modality]
                scaled_data = scaler.transform(data)
            
            self.omic_data[modality] = scaled_data


    def _process_stages(self, stages):
        """Convert stage strings to numeric labels."""
        stage_mapping = {'I': 0, 'II': 1, 'III': 2, 'IV': 3}
        
        def extract_main_stage(stage_str):
            # Handle special cases
            if 'I/II' in stage_str:
                self.logger.warning(f"Found ambiguous stage: {stage_str}, treating as stage II")
                return 'II'
            if 'II/III' in stage_str:
                self.logger.warning(f"Found ambiguous stage: {stage_str}, treating as stage III")
                return 'III'
            
            # Remove 'Stage ' prefix if it exists
            if 'Stage ' in stage_str:
                stage_str = stage_str.replace('Stage ', '')
                
            # Check stages in reverse order to avoid false matches
            if stage_str.startswith('IV'):
                return 'IV'
            elif stage_str.startswith('III'):
                return 'III'
            elif stage_str.startswith('II'):
                return 'II'
            elif stage_str.startswith('I'):
                return 'I'
                
            return stage_str

        processed_stages = []
        skipped_stages = 0
        
        # Print all unique stages before processing
        unique_stages = set(stages)
        self.logger.info(f"Unique stages in data: {unique_stages}")
        
        for stage in stages:
            try:
                main_stage = extract_main_stage(stage)
                if main_stage not in stage_mapping:
                    self.logger.warning(f"Skipping unknown stage: {stage} (extracted: {main_stage})")
                    skipped_stages += 1
                    continue
                processed_stages.append(stage_mapping[main_stage])
            except Exception as e:
                self.logger.warning(f"Error processing stage '{stage}': {str(e)}")
                skipped_stages += 1
                continue
        
        # Log stage distribution
        unique, counts = np.unique(processed_stages, return_counts=True)
        dist_dict = dict(zip([f"Stage {list(stage_mapping.keys())[v]}" for v in unique], counts))
        self.logger.info(f"Stage distribution after processing: {dist_dict}")
        if skipped_stages > 0:
            self.logger.warning(f"Skipped {skipped_stages} invalid stages")
        
        return processed_stages

    def __len__(self):
        return len(self.clinical_data)

    def __getitem__(self, idx):
        # Get omic data for each modality
        sample = {
            modality: torch.FloatTensor(data[idx])
            for modality, data in self.omic_data.items()
        }
        
        survival_time = torch.FloatTensor([self.survival_times[idx]])
        stage_label = torch.LongTensor([self.stage_labels[idx]])
        
        return sample, survival_time, stage_label

    def get_scalers(self):
        """Return fitted scalers for use with validation/test data."""
        scalers = self.scalers.copy()
        scalers['survival'] = self.survival_scaler
        return scalers