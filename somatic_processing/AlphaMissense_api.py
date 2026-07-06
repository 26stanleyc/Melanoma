import pandas as pd
import aiohttp
import asyncio
import time
import json
import os
from tqdm import tqdm

def save_checkpoint(completed_variants, pathogenic_variants, checkpoint_file):
   """Save progress to checkpoint file"""
   checkpoint = {
       'completed_variants': completed_variants,
       'pathogenic_variants': pathogenic_variants
   }
   with open(checkpoint_file, 'w') as f:
       json.dump(checkpoint, f)

def load_checkpoint(checkpoint_file):
   """Load progress from checkpoint file"""
   if os.path.exists(checkpoint_file):
       with open(checkpoint_file, 'r') as f:
           return json.load(f)
   return {'completed_variants': [], 'pathogenic_variants': []}

def parse_aa_string(aa_string):
   """Parse strings like '2:F,Y' into list of amino acids"""
   if not aa_string or aa_string == '"':  # Handle empty cases
       return []
   try:
       count, aas = aa_string.split(':')
       return aas.split(',')
   except ValueError:
       return []

async def fetch_data_with_retry(session, protein, position, original_aa, semaphore, max_retries=10):
   """Fetch data with retries and amino acid validation"""
   url = "https://alphamissense.hegelab.org/hotspotapi"
   params = {"uid": protein, "resi": position}
   full_url = f"{url}?uid={protein}&resi={position}"
   
   for attempt in range(max_retries):
       async with semaphore:
           try:
               async with session.get(url, params=params, timeout=30) as response:
                   if response.status == 404:
                       print(f"404 Error - Resource not found for {protein}:{position}")
                       error_info = {
                           'protein': protein,
                           'position': position,
                           'original_aa': original_aa,
                           'error': '404 - Resource not found',
                           'url': full_url,
                           'attempts': 1
                       }
                       return None, error_info
                       
                   if response.status == 200:
                       data = await response.json()
                       
                       # Check if the position exists and reference AA matches
                       if 'aa' not in data:
                           print(f"Position {position} not found in {protein}")
                           return None, None
                           
                       if data['aa'] != original_aa:
                           print(f"Reference AA mismatch for {protein}:{position} - Expected {original_aa}, got {data['aa']}")
                           return None, None
                           
                       data['pathogenic_aas'] = parse_aa_string(data.get('pathogenic', ''))
                       return data, None
                       
                   print(f"Attempt {attempt+1} failed with status {response.status} for URL: {full_url}")
               await asyncio.sleep(0.5)
           except Exception as e:
               print(f"Attempt {attempt+1} failed for URL: {full_url}")
               print(f"Error: {str(e)}")
               if attempt == max_retries - 1:
                   error_info = {
                       'protein': protein,
                       'position': position,
                       'original_aa': original_aa,
                       'error': str(e),
                       'url': full_url,
                       'attempts': attempt + 1
                   }
                   return None, error_info
               await asyncio.sleep(0.5)
               continue
   return None, None

def save_errors(errors, error_file):
   """Save failed queries to file"""
   error_df = pd.DataFrame(errors)
   if os.path.exists(error_file):
       existing_errors = pd.read_csv(error_file)
       error_df = pd.concat([existing_errors, error_df], ignore_index=True)
   error_df.to_csv(error_file, index=False)

async def process_variants(input_file, output_file, error_file='failed_queries.csv', 
                        checkpoint_file='checkpoint.json', max_concurrent=100):
   """Process variants with retry logic and error logging"""
   start_time = time.time()
   # Load checkpoint
   checkpoint = load_checkpoint(checkpoint_file)
   completed_variants = set(checkpoint['completed_variants'])
   pathogenic_variants = checkpoint['pathogenic_variants']
   errors = []
   
   # Read and filter data
   df = pd.read_csv(input_file, sep=',')
   remaining_df = df[~df.index.astype(str).isin(completed_variants)]
   
   print(f"Starting processing of {len(remaining_df)} variants...")
   
   # Semaphore for rate limiting
   semaphore = asyncio.Semaphore(max_concurrent)
   
   try:
       async with aiohttp.ClientSession() as session:
           # Process in smaller batches
           batch_size = 100
           for start_idx in range(0, len(remaining_df), batch_size):
               batch_df = remaining_df.iloc[start_idx:start_idx + batch_size]
               print(f"Processing batch starting at index {start_idx}")
               
               # Create tasks for current batch
               tasks = []
               for idx, row in batch_df.iterrows():
                   print(f"Creating task for {row['Protein']}:{row['Position']} Original:{row['Original_AA']}")
                   task = asyncio.create_task(
                       fetch_data_with_retry(session, row['Protein'], row['Position'], row['Original_AA'], semaphore)
                   )
                   tasks.append((idx, row, task))
               
               # Process batch results
               with tqdm(total=len(tasks), desc=f"Batch {start_idx//batch_size + 1}") as pbar:
                   for i, (idx, row, task) in enumerate(tasks):
                       try:
                           data, error_info = await task
                           if error_info:
                               errors.append(error_info)
                               print(f"Failed: {row['Protein']} {row['Position']}")
                           elif data and 'pathogenic_aas' in data:
                               if row['New_AA'] in data['pathogenic_aas']:
                                   pathogenic_variants.append(row.to_dict())
                                   print(f"Found pathogenic: {row['Protein']} {row['Position']}{row['New_AA']}")
                           
                           completed_variants.add(str(idx))
                           pbar.update(1)
                           
                           # Save checkpoint and errors periodically
                           if (i + 1) % 50 == 0:
                               save_checkpoint(list(completed_variants), pathogenic_variants, checkpoint_file)
                               if errors:
                                   save_errors(errors, error_file)
                           
                       except Exception as e:
                           print(f"Error processing {row['Protein']} {row['Position']}: {e}")
               
               # Save checkpoint after each batch
               save_checkpoint(list(completed_variants), pathogenic_variants, checkpoint_file)
               if errors:
                   save_errors(errors, error_file)
   
   except KeyboardInterrupt:
       print("\nInterrupt detected. Saving checkpoint and errors...")
       save_checkpoint(list(completed_variants), pathogenic_variants, checkpoint_file)
       if errors:
           save_errors(errors, error_file)
       print("Checkpoint and errors saved. You can resume later.")
       return
   
   # Save final results
   if pathogenic_variants:
       pathogenic_df = pd.DataFrame(pathogenic_variants)
       pathogenic_df.to_csv(output_file, index=False)
       print(f"\nSaved {len(pathogenic_variants)} pathogenic variants to {output_file}")
   
   # Save final errors
   if errors:
       save_errors(errors, error_file)
       print(f"\nSaved {len(errors)} failed queries to {error_file}")
   
   # Clean up checkpoint after successful completion
   if os.path.exists(checkpoint_file):
       os.remove(checkpoint_file)
   
   # Print final statistics
   elapsed = time.time() - start_time
   final_rate = len(remaining_df) / elapsed
   print(f"\nFinal Statistics:")
   print(f"Total variants processed: {len(df)}")
   print(f"Failed queries: {len(errors)}")
   print(f"Average processing rate: {final_rate:.2f} variants/second")

if __name__ == "__main__":
   input_file = "/Users/stanleychen/git/Melanoma/data/somatic_mutation_processed_v2/all_missense_mutations.csv"
   output_file = "pathogenic_variants.csv"
   error_file = "failed_queries.csv"
   checkpoint_file = "alphamissense_checkpoint.json"
   
   asyncio.run(process_variants(input_file, output_file, error_file, checkpoint_file))