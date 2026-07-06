import pandas as pd

def compare_clinical_patients(file1_path, file2_path):
    """
    Compare patient IDs between two clinical data files.
    
    Args:
        file1_path (str): Path to first clinical data file
        file2_path (str): Path to second clinical data file
    """
    try:
        # Read both files
        df1 = pd.read_csv(file1_path)
        df2 = pd.read_csv(file2_path)
        
        # Get patient IDs from first column
        patients1 = set(df1['patient_barcode'])
        patients2 = set(df2['patient_barcode'])
        
        # Find overlapping IDs
        overlapping = patients1.intersection(patients2)
        
        # Print results
        print("\nResults:")
        print(f"Patients in first file: {len(patients1)}")
        print(f"Patients in second file: {len(patients2)}")
        print(f"Number of overlapping patients: {len(overlapping)}")
        
        # Print some example overlapping IDs
        print("\nFirst 5 overlapping patient IDs:")
        for id in list(overlapping)[:5]:
            print(f"- {id}")
            
        # Print unique to each file
        unique_to_first = patients1 - patients2
        unique_to_second = patients2 - patients1
        
        print(f"\nUnique to first file: {len(unique_to_first)} patients")
        print("Examples:")
        for id in list(unique_to_first)[:3]:
            print(f"- {id}")
            
        print(f"\nUnique to second file: {len(unique_to_second)} patients")
        print("Examples:")
        for id in list(unique_to_second)[:3]:
            print(f"- {id}")
        
        # # Save results to file
        # with open('patient_overlap_results.txt', 'w') as f:
        #     f.write("Patient ID Overlap Analysis\n\n")
        #     f.write(f"Patients in first file: {len(patients1)}\n")
        #     f.write(f"Patients in second file: {len(patients2)}\n")
        #     f.write(f"Overlapping patients: {len(overlapping)}\n\n")
        #     f.write("Overlapping IDs:\n")
        #     for id in sorted(overlapping):
        #         f.write(f"{id}\n")
        
        return {
            'overlapping': overlapping,
            'unique_to_first': unique_to_first,
            'unique_to_second': unique_to_second
        }
        
    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}")
    except Exception as e:
        print(f"Error: An unexpected error occurred - {e}")

# File paths - replace these with your actual file paths
FILE1_PATH = "/Users/stanleychen/git/Melanoma/final_data/filtered_patient_data_all_v2.csv"
FILE2_PATH = "/Users/stanleychen/git/Melanoma/final_data/final_label_v2.csv"
compare_clinical_patients(FILE1_PATH, FILE2_PATH)

