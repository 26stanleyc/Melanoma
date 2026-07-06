import pandas as pd
df = pd.read_csv('/Users/stanleychen/git/Melanoma/combined_RNA_data.csv', sep=',', nrows=1)  # Read only first row
print(f"Number of columns: {len(df.columns)}")