def count_tumor_stages_by_split(file_path):
    # Initialize counters for each stage and split
    train_counts = {0: 0, 1: 0, 2: 0, 3: 0}
    validation_counts = {0: 0, 1: 0, 2: 0, 3: 0}
    
    # Read the file
    with open(file_path, 'r') as file:
        # Skip header
        header = file.readline()
        
        # Process each line
        for line in file:
            # Split the line by comma
            fields = line.strip().split(',')
            
            try:
                # Get tumor stage (last field)
                stage = int(float(fields[-1]))
                
                # Get split (train/validation)
                split = fields[4]
                
                # Count based on split
                if split == 'train':
                    if stage in train_counts:
                        train_counts[stage] += 1
                elif split == 'validation':
                    if stage in validation_counts:
                        validation_counts[stage] += 1
                
            except (ValueError, IndexError):
                continue
    
    # Print results
    print("\nTraining set samples:")
    for stage in sorted(train_counts.keys()):
        print(f"Stage {stage}: {train_counts[stage]} samples")
    
    print("\nValidation set samples:")
    for stage in sorted(validation_counts.keys()):
        print(f"Stage {stage}: {validation_counts[stage]} samples")
    
    return train_counts, validation_counts

# Example usage:
train_counts, validation_counts = count_tumor_stages_by_split('/Users/stanleychen/git/Melanoma/final_data/final_label_v2.csv')