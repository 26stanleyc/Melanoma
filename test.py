import pandas as pd

# Read the CSV file
df = pd.read_csv('/Users/stanleychen/git/Melanoma/multiomics-open-research/data/transformed_data.csv')  # Replace with your file path

# Get dimensions
rows, cols = df.shape

# Print basic dimension info
print(f"Dimensions of the CSV file:")
print(f"Number of rows: {rows}")
print(f"Number of columns: {cols}")

# Print additional info
print("\nDetailed information:")
print(f"Memory usage: {df.memory_usage().sum() / 1024**2:.2f} MB")
print("\nColumn names:")
print(df.columns.tolist())

# Print first few rows to verify structure
print("\nFirst few rows:")
print(df.head())