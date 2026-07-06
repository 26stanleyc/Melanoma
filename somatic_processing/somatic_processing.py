import pandas as pd
import re
import os
from pathlib import Path

def parse_maf_file(file_path):
    """
    Parse MAF file with proper handling of headers and comments.
    """
    print(f"Reading MAF file: {file_path}")
    
    # Read MAF file while skipping comments
    df = pd.read_csv(file_path, sep='\t', comment='#')
    print(f"Found {len(df)} total variants")
    
    # Filter for missense mutations
    missense_df = df[df['Variant_Classification'] == 'Missense_Mutation'].copy()
    print(f"Found {len(missense_df)} missense mutations")
    
    result = []
    for _, row in missense_df.iterrows():
        protein = row['Hugo_Symbol']
        
        if pd.notna(row['HGVSp_Short']):
            match = re.match(r'p\.([A-Z])(\d+)([A-Z])', row['HGVSp_Short'])
            if match:
                original_aa, position, new_aa = match.groups()
                
                result.append({
                    'Protein': protein,
                    'Original_AA': original_aa,
                    'Position': int(position),
                    'New_AA': new_aa,
                    'Chromosome': row['Chromosome'],
                    'Sample_ID': row['Tumor_Sample_Barcode'],
                    'Transcript': row['Transcript_ID'],
                    'SIFT': row['SIFT'],
                    'PolyPhen': row['PolyPhen']
                })
    
    result_df = pd.DataFrame(result)
    print(f"Processed {len(result_df)} mutations with valid amino acid changes")
    
    if not result_df.empty:
        result_df['Mutation'] = result_df.apply(
            lambda x: f"{x['Protein']} {x['Original_AA']}{x['Position']}{x['New_AA']}", 
            axis=1
        )
        
        columns = ['Mutation', 'Protein', 'Original_AA', 'Position', 'New_AA', 
                  'Chromosome', 'Sample_ID', 'Transcript', 'SIFT', 'PolyPhen']
        
        return result_df[columns]
    
    return pd.DataFrame()

def process_mutation_directory(input_dir, output_dir):
    """
    Process all MAF files in the input directory and save results to output directory.
    """
    # Convert to absolute paths
    input_dir = os.path.expanduser(input_dir)
    output_dir = os.path.expanduser(output_dir)
    
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize a list to store all mutations
    all_mutations = []
    
    # Process only MAF files
    input_path = Path(input_dir)
    
    files_found = False
    for file_path in input_path.glob('*.maf'):
        files_found = True
        try:
            print(f"\nProcessing {file_path.name}...")
            mutations_df = parse_maf_file(file_path)
            
            if not mutations_df.empty:
                mutations_df['Source_File'] = file_path.name
                output_file = Path(output_dir) / f"processed_{file_path.stem}.csv"
                mutations_df.to_csv(output_file, index=False)
                print(f"Saved results to {output_file}")
                all_mutations.append(mutations_df)
            else:
                print(f"No valid mutations found in {file_path.name}")
        
        except Exception as e:
            print(f"Error processing {file_path.name}: {str(e)}")
            print("Full error details:", e.__class__.__name__)
            import traceback
            traceback.print_exc()
    
    if not files_found:
        print("\nNo MAF files found in directory")
        print("Directory contents:")
        for item in input_path.iterdir():
            print(f"- {item.name}")
    
    if all_mutations:
        combined_df = pd.concat(all_mutations, ignore_index=True)
        combined_output = Path(output_dir) / "all_missense_mutations.csv"
        combined_df.to_csv(combined_output, index=False)
        print(f"\nSaved combined results to {combined_output}")
        print(f"Total missense mutations found: {len(combined_df)}")
    else:
        print("\nNo mutations were processed successfully")

if __name__ == "__main__":
    input_directory = "/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed"
    output_directory = "/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed_v2"
    
    print("Starting mutation processing script...")
    process_mutation_directory(input_directory, output_directory)
    print("Script completed")