import pandas as pd
import os
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
                           usecols=['#CHROM', 'POS', 'REF', 'ALT', 'am_pathogenicity', 'am_class'],
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
        return alpha_df.iloc[idx]['am_pathogenicity'], alpha_df.iloc[idx]['am_class']
    return None, None

def merge_alphamissense_scores(mutations_file, alphamissense_file, output_file, checkpoint_dir='checkpoints'):
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
        found_scores = checkpoint['found_scores']
        print(f"Resuming from mutation {start_idx}")
    else:
        print("\nReading mutations file...")
        mutations_df = pd.read_csv(mutations_file, sep=',', low_memory=False)
        start_idx = 0
        found_scores = 0
        
        # Add AlphaMissense score columns
        mutations_df['AlphaMissense_Score'] = pd.NA
        mutations_df['AlphaMissense_Class'] = pd.NA
    
    print(f"Processing {len(mutations_df) - start_idx} remaining mutations...")
    
    checkpoint_interval = 1000  # Save checkpoint every 1000 mutations
    
    for idx in range(start_idx, len(mutations_df)):
        if idx % 100 == 0:
            print(f"Processing mutation {idx}/{len(mutations_df)}")
        
        row = mutations_df.iloc[idx]
        
        score, pred_class = binary_search_alphamissense(
            alpha_df, 
            row['chromosome'],
            row['position'],
            row['reference_allele'],
            row['variant_allele']
        )
        
        if score is not None:
            mutations_df.loc[idx, 'AlphaMissense_Score'] = score
            mutations_df.loc[idx, 'AlphaMissense_Class'] = pred_class
            found_scores += 1
            
            if found_scores <= 5:
                print(f"\nMatch found:")
                print(f"Location: {row['chromosome']}:{row['position']} {row['reference_allele']}>{row['variant_allele']}")
                print(f"Score: {score}")
                print(f"Class: {pred_class}")
        
        # Save checkpoint periodically
        if (idx + 1) % checkpoint_interval == 0:
            checkpoint_state = {
                'mutations_df': mutations_df,
                'current_idx': idx + 1,
                'found_scores': found_scores
            }
            save_checkpoint(checkpoint_state, progress_file)
    
    # Save final results
    print(f"\nSaving results to {output_file}")
    mutations_df.to_csv(output_file, index=False)
    
    # Clean up checkpoint files
    if os.path.exists(progress_file):
        os.remove(progress_file)
    
    print(f"\nSummary:")
    print(f"Total mutations: {len(mutations_df)}")
    print(f"Found AlphaMissense scores: {found_scores}")
    print(f"Coverage: {found_scores/len(mutations_df)*100:.1f}%")

if __name__ == "__main__":
    merge_alphamissense_scores(
        mutations_file="/Users/stanleychen/git/Melanoma/final_data/missense_mutations.csv",
        alphamissense_file="/Users/stanleychen/git/Melanoma/somatic_processing/AlphaMissense_hg38.tsv",
        output_file="mutations_with_alphamissense.csv",
        checkpoint_dir="checkpoints"
    )