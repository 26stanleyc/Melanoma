import pandas as pd
import os
from pathlib import Path
import numpy as np
from bisect import bisect_left
import pickle

def save_checkpoint(state, filename):
    """Save current state to checkpoint file"""
    print(f"Saving checkpoint to {filename}...")
    with open(filename, 'wb') as f:
        pickle.dump(state, f)

def load_checkpoint(filename):
    """Load state from checkpoint file"""
    if os.path.exists(filename):
        print(f"Loading checkpoint from {filename}...")
        with open(filename, 'rb') as f:
            return pickle.load(f)
    return None

def create_alphamissense_index(alphamissense_file, checkpoint_dir):
    """Create sorted index of AlphaMissense file for fast lookups"""
    checkpoint_file = os.path.join(checkpoint_dir, 'alphamissense_index.pkl')
    
    # Try to load from checkpoint
    alpha_df = load_checkpoint(checkpoint_file)
    if alpha_df is not None:
        print(f"Loaded indexed AlphaMissense with {len(alpha_df):,} entries")
        return alpha_df
    
    print("Creating AlphaMissense index...")
    chunks = []
    total_rows = 0
    
    for chunk in pd.read_csv(alphamissense_file, 
                           sep='\t', 
                           skiprows=3, 
                           usecols=['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity'],
                           low_memory=False,
                           chunksize=100000):
        chunk['lookup_key'] = chunk.apply(lambda x: f"{x['#CHROM']}_{x['POS']}_{x['REF']}_{x['ALT']}", axis=1)
        chunks.append(chunk)
        total_rows += len(chunk)
        print(f"Processed {total_rows:,} rows...")
    
    alpha_df = pd.concat(chunks)
    print("Sorting index...")
    alpha_df.sort_values('lookup_key', inplace=True)
    
    # Save index checkpoint
    save_checkpoint(alpha_df, checkpoint_file)
    
    return alpha_df

def binary_search_alphamissense(alpha_df, chrom, pos, ref, alt):
    """Binary search for matching mutation in sorted AlphaMissense data"""
    search_key = f"{chrom}_{pos}_{ref}_{alt}"
    idx = alpha_df['lookup_key'].searchsorted(search_key)
    
    if idx < len(alpha_df) and alpha_df.iloc[idx]['lookup_key'] == search_key:
        return alpha_df.iloc[idx]['am_pathogenicity']
    return None

def read_processed_maf_file(filepath):
    """Read MAF file from processed directory"""
    try:
        needed_columns = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 
                         'Reference_Allele', 'Tumor_Seq_Allele2', 
                         'Transcript_ID', 'Variant_Classification']
        
        df = pd.read_csv(filepath, 
                        sep='\t', 
                        comment='#',
                        low_memory=False,
                        usecols=needed_columns)
        return df
    except Exception as e:
        print(f"Error reading processed MAF file {filepath}: {e}")
        return None

def find_coordinates_in_processed_files(mutation_row, processed_dir, coordinate_cache=None):
    """Find genomic coordinates for a mutation by searching processed files"""
    # Check cache first
    if coordinate_cache is not None:
        cache_key = f"{mutation_row['Mutation']}_{mutation_row['Transcript']}"
        if cache_key in coordinate_cache:
            return coordinate_cache[cache_key]
    
    gene_name = mutation_row['Mutation'].split()[0]
    transcript = mutation_row['Transcript']
    
    for filename in os.listdir(processed_dir):
        if not filename.endswith('.maf'):
            continue
            
        processed_df = read_processed_maf_file(os.path.join(processed_dir, filename))
        if processed_df is None:
            continue
            
        try:
            matches = processed_df[
                (processed_df['Hugo_Symbol'] == gene_name) &
                (processed_df['Transcript_ID'] == transcript)
            ]
            
            if len(matches) > 0:
                match = matches.iloc[0]
                result = {
                    'chrom': match['Chromosome'],
                    'pos': match['Start_Position'],
                    'ref': match['Reference_Allele'],
                    'alt': match['Tumor_Seq_Allele2'],
                    'gene': gene_name,
                    'variant_class': match['Variant_Classification'],
                    'source_file': filename
                }
                
                # Add to cache
                if coordinate_cache is not None:
                    coordinate_cache[cache_key] = result
                    
                return result
                
        except Exception as e:
            print(f"Error processing matches in file {filename}: {e}")
            continue
    
    return None

def merge_alphamissense_scores(mutations_file, processed_dir, alphamissense_file, output_file, checkpoint_dir='checkpoints'):
    """Add AlphaMissense scores to mutations using indexed lookups with checkpoints"""
    # Create checkpoint directory if it doesn't exist
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    # Load or create AlphaMissense index
    alpha_df = create_alphamissense_index(alphamissense_file, checkpoint_dir)
    print(f"Using index with {len(alpha_df):,} entries")
    
    # Try to load progress checkpoint
    progress_file = os.path.join(checkpoint_dir, 'merge_progress.pkl')
    checkpoint = load_checkpoint(progress_file)
    
    if checkpoint is not None:
        mutations_df = checkpoint['mutations_df']
        start_idx = checkpoint['current_idx']
        found_coords = checkpoint['found_coords']
        found_scores = checkpoint['found_scores']
        coordinate_cache = checkpoint['coordinate_cache']
        print(f"Resuming from mutation {start_idx}")
    else:
        print("\nReading mutations file...")
        # Read original mutations file and store its columns
        mutations_df = pd.read_csv(mutations_file, sep=',', low_memory=False)
        original_columns = mutations_df.columns.tolist()
        start_idx = 0
        found_coords = 0
        found_scores = 0
        coordinate_cache = {}
        
        # Add only AlphaMissense score column
        mutations_df['AlphaMissense_Score'] = pd.NA
    
    print(f"Processing {len(mutations_df) - start_idx} remaining mutations...")
    
    checkpoint_interval = 1000  # Save checkpoint every 1000 mutations
    
    for idx in range(start_idx, len(mutations_df)):
        if idx % 100 == 0:
            print(f"Processing mutation {idx}/{len(mutations_df)}")
        
        row = mutations_df.iloc[idx]
        coords = find_coordinates_in_processed_files(row, processed_dir, coordinate_cache)
        
        if coords:
            found_coords += 1
            
            score = binary_search_alphamissense(
                alpha_df, 
                coords['chrom'], 
                coords['pos'], 
                coords['ref'], 
                coords['alt']
            )
            
            if score is not None:
                mutations_df.loc[idx, 'AlphaMissense_Score'] = score
                found_scores += 1
                
                if found_scores <= 5:
                    print(f"\nMatch found:")
                    print(f"Original: {row['Mutation']}")
                    print(f"Location: {coords['chrom']}:{coords['pos']} {coords['ref']}>{coords['alt']}")
                    print(f"Score: {score}")
        
        # Save checkpoint periodically
        if (idx + 1) % checkpoint_interval == 0:
            checkpoint_state = {
                'mutations_df': mutations_df,
                'current_idx': idx + 1,
                'found_coords': found_coords,
                'found_scores': found_scores,
                'coordinate_cache': coordinate_cache
            }
            save_checkpoint(checkpoint_state, progress_file)
    
    # Ensure only original columns plus AlphaMissense_Score are in the output
    final_columns = original_columns + ['AlphaMissense_Score']
    mutations_df = mutations_df[final_columns]
    
    # Save final results
    print(f"\nSaving results to {output_file}")
    mutations_df.to_csv(output_file, index=False)
    
    # Clean up checkpoint files
    if os.path.exists(progress_file):
        os.remove(progress_file)
    
    print(f"\nSummary:")
    print(f"Total mutations: {len(mutations_df)}")
    print(f"Found coordinates: {found_coords}")
    print(f"Found AlphaMissense scores: {found_scores}")
    print(f"Coverage: {found_scores/len(mutations_df)*100:.1f}%")

if __name__ == "__main__":
    merge_alphamissense_scores(
        mutations_file="/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/maf/pathogenic_variants.csv",
        processed_dir="/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed",
        alphamissense_file="/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed_v3/AlphaMissense_hg38.tsv",
        output_file="mutations_with_alphamissense.csv",
        checkpoint_dir="/Users/stanleychen/git/Melanoma/final_data/not_normalized_data/maf/checkpoints_for_am"
    )