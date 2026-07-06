import pandas as pd
import os
from pathlib import Path
import numpy as np

def read_maf_alleles(filepath, gene, pos):
    """Read ref/alt alleles from MAF file for specific position"""
    try:
        df = pd.read_csv(filepath, 
                        sep='\t',
                        comment='#',
                        usecols=['Hugo_Symbol', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'],
                        low_memory=False)
        
        matches = df[
            (df['Hugo_Symbol'] == gene) &
            (df['Start_Position'] == pos)
        ]
        
        if len(matches) > 0:
            match = matches.iloc[0]
            return match['Reference_Allele'], match['Tumor_Seq_Allele2']
    except Exception as e:
        print(f"Error reading MAF file {filepath}: {e}")
    return None, None

def create_alphamissense_index(alphamissense_file):
    """Create sorted index of AlphaMissense file for fast lookups"""
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
    
    return alpha_df

def binary_search_alphamissense(alpha_df, chrom, pos, ref, alt):
    """Binary search for matching mutation in sorted AlphaMissense data"""
    search_key = f"{chrom}_{pos}_{ref}_{alt}"
    idx = alpha_df['lookup_key'].searchsorted(search_key)
    
    if idx < len(alpha_df) and alpha_df.iloc[idx]['lookup_key'] == search_key:
        return alpha_df.iloc[idx]['am_pathogenicity']
    return None

def add_alphamissense_scores(mutations_file, processed_dir, alphamissense_file, output_file):
    """Add AlphaMissense scores using pre-computed genomic coordinates"""
    
    # Create AlphaMissense index
    alpha_df = create_alphamissense_index(alphamissense_file)
    print(f"Created index with {len(alpha_df):,} entries")
    
    # Read mutations file
    print("\nReading mutations with coordinates...")
    mutations_df = pd.read_csv(mutations_file, low_memory=False)
    print(f"Loaded {len(mutations_df)} mutations")
    
    # Cache for MAF file contents
    maf_cache = {}
    
    # Process each mutation using genomic coordinates
    print("\nLooking up AlphaMissense scores...")
    found_scores = 0
    
    # Add new column for scores if it doesn't exist
    if 'AlphaMissense_Score' not in mutations_df.columns:
        mutations_df['AlphaMissense_Score'] = pd.NA
    
    for idx, row in mutations_df.iterrows():
        if idx % 1000 == 0:
            print(f"Processed {idx}/{len(mutations_df)} mutations...")
            
        try:
            if pd.notna(row['Genomic_Position']) and pd.notna(row['Source_MAF_File']):
                # Parse genomic coordinate
                chrom, pos = row['Genomic_Position'].split(':')
                pos = int(pos)
                gene = row['Mutation'].split()[0]
                
                # Get ref/alt alleles from source MAF file
                maf_path = os.path.join(processed_dir, row['Source_MAF_File'])
                
                if maf_path not in maf_cache:
                    # Read and cache MAF file contents
                    ref, alt = read_maf_alleles(maf_path, gene, pos)
                    if ref is not None and alt is not None:
                        maf_cache[maf_path] = (ref, alt)
                else:
                    ref, alt = maf_cache[maf_path]
                
                if ref is not None and alt is not None:
                    score = binary_search_alphamissense(alpha_df, chrom, pos, ref, alt)
                    
                    if score is not None:
                        mutations_df.loc[idx, 'AlphaMissense_Score'] = score
                        found_scores += 1
                        
                        if found_scores <= 5:
                            print(f"\nMatch found:")
                            print(f"Gene: {gene}")
                            print(f"Location: {row['Genomic_Position']}")
                            print(f"Ref/Alt: {ref}/{alt}")
                            print(f"Score: {score}")
                            print(f"Source: {row['Source_MAF_File']}")
        except Exception as e:
            print(f"Error processing mutation {idx}: {e}")
            continue
    
    # Save results
    print(f"\nSaving results to {output_file}")
    mutations_df.to_csv(output_file, index=False)
    
    # Print summary
    print(f"\nSummary:")
    print(f"Total mutations: {len(mutations_df)}")
    print(f"Found AlphaMissense scores: {found_scores}")
    print(f"Coverage: {found_scores/len(mutations_df)*100:.1f}%")

if __name__ == "__main__":
    add_alphamissense_scores(
        mutations_file="/Users/stanleychen/git/Melanoma/mutations_with_alphamissense.csv",
        processed_dir="/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed",
        alphamissense_file="/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed_v3/AlphaMissense_hg38.tsv",
        output_file="mutations_with_scores_final.csv"
    )