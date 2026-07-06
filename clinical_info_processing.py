import csv
from collections import defaultdict

def process_clinical_drug_data(input_file, output_file):
    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        data = list(reader)

    # Count [Not Available] and [Not Applicable] in each column
    not_available_count = defaultdict(int)
    total_rows = len(data)
    for row in data:
        for i, value in enumerate(row):
            if value in ["[Not Available]", "[Not Applicable]"]:
                not_available_count[i] += 1

    # Identify columns to keep
    columns_to_keep = [i for i in range(len(headers)) if not_available_count[i] <= total_rows / 2]

    # Ensure important identifier columns are always kept
    important_columns = ['bcr_patient_uuid', 'bcr_patient_barcode', 'bcr_drug_barcode', 'bcr_drug_uuid']
    for col in important_columns:
        if col in headers and headers.index(col) not in columns_to_keep:
            columns_to_keep.append(headers.index(col))
    columns_to_keep.sort()

    # Create new headers and data
    new_headers = [headers[i] for i in columns_to_keep]
    new_data = [[row[i] for i in columns_to_keep] for row in data]

    # Write to new file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(new_headers)
        writer.writerows(new_data)

    # Identify removed columns
    removed_columns = [headers[i] for i in range(len(headers)) if i not in columns_to_keep]

    return removed_columns, new_headers

# Process the data
input_file = '/Users/stanleychen/git/Melanoma/clinical_data/clinical_tme_info/nationwidechildrens.org_clinical_nte_skcm.txt'  # Assume we saved the input data to this file
output_file = '/Users/stanleychen/git/Melanoma/clinical_data_processed/clinical_tme_processed.tsv'
removed_columns, kept_columns = process_clinical_drug_data(input_file, output_file)

print("Columns removed due to more than 50% [Not Available] or [Not Applicable]:")
for column in removed_columns:
    print(f"- {column}")

print("\nColumns kept in the processed file:")
for column in kept_columns:
    print(f"- {column}")