import os
import pandas as pd
from typing import Dict, List

class MAFParser:
    def __init__(self):
        self.mutation_info = []
        
    def process_directory(self, root_dir: str) -> List[Dict]:
        """
        Walk through directory and process all .maf files (excluding .gz files)
        
        Args:
            root_dir (str): Root directory path to start searching
            
        Returns:
            List[Dict]: List of mutation information dictionaries
        """
        for root, _, files in os.walk(root_dir):
            for file in files:
                if file.endswith('.maf') and not file.endswith('.gz'):
                    file_path = os.path.join(root, file)
                    self.process_maf_file(file_path)
        return self.mutation_info

    def process_maf_file(self, file_path: str) -> None:
        """
        Process individual MAF file and extract mutation information
        
        Args:
            file_path (str): Path to MAF file
        """
        try:
            # Read MAF file using pandas
            df = pd.read_csv(file_path, sep='\t', comment='#', low_memory=False)
            
            # Filter for only missense mutations
            missense_df = df[df['Variant_Classification'] == 'Missense_Mutation']
            
            for _, row in missense_df.iterrows():
                mutation_dict = {
                    'gene': row['Hugo_Symbol'],
                    'chromosome': row['Chromosome'],
                    'position': row['Start_Position'],
                    'reference_allele': row['Reference_Allele'],
                    'variant_allele': row['Tumor_Seq_Allele2'],
                    'protein_change': row['HGVSp_Short'],
                    'nucleotide_change': row['HGVSc'],
                    'variant_type': row['Variant_Classification'],
                    'amino_acid_change': self._parse_amino_acid_change(row['HGVSp_Short']),
                    'sift': row['SIFT'],
                    'polyphen': row['PolyPhen']
                }
                self.mutation_info.append(mutation_dict)
                
        except Exception as e:
            print(f"Error processing file {file_path}: {str(e)}")
    
    def _parse_amino_acid_change(self, hgvsp: str) -> Dict[str, str]:
        """
        Parse the amino acid change from HGVSp notation
        
        Args:
            hgvsp (str): HGVSp notation string (e.g., "p.R645W")
            
        Returns:
            Dict[str, str]: Dictionary containing original and new amino acids and position
        """
        if not hgvsp or pd.isna(hgvsp):
            return {'original': '', 'new': '', 'position': ''}
            
        # Remove the 'p.' prefix if present
        if hgvsp.startswith('p.'):
            hgvsp = hgvsp[2:]
            
        try:
            original = hgvsp[0]
            position = ''.join(filter(str.isdigit, hgvsp))
            new = hgvsp[-1]
            
            return {
                'original': original,
                'new': new,
                'position': position
            }
        except:
            return {'original': '', 'new': '', 'position': ''}

def save_mutations_to_csv(mutations: List[Dict], output_file: str) -> None:
    """
    Save mutation information to a CSV file
    
    Args:
        mutations (List[Dict]): List of mutation dictionaries
        output_file (str): Path to output CSV file
    """
    # Convert the nested amino acid change dictionary into separate columns
    formatted_mutations = []
    for mut in mutations:
        mut_copy = mut.copy()
        amino_acids = mut_copy.pop('amino_acid_change')
        mut_copy.update({
            'original_aa': amino_acids['original'],
            'new_aa': amino_acids['new'],
            'aa_position': amino_acids['position']
        })
        formatted_mutations.append(mut_copy)
    
    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(formatted_mutations)
    df.to_csv(output_file, index=False)
    print(f"Saved {len(mutations)} mutations to {output_file}")

def main():
    # Example usage
    parser = MAFParser()
    
    # Process directory and get mutation information
    mutations = parser.process_directory('/Users/stanleychen/git/Melanoma/data/maf')
    
    # Save to CSV file
    output_file = "missense_mutations.csv"
    save_mutations_to_csv(mutations, output_file)

if __name__ == "__main__":
    main()