import pandas as pd

# Read the file
# Assuming the second file is a CSV or similar, first let's read it properly
df = pd.read_csv('/Users/stanleychen/git/Melanoma/multiomics-open-research/data/transformed_data.csv')  # Replace with your actual file path

# Get all the ENSG IDs (column names) and join them into one line
ensg_line = ','.join(df.columns)

# Get the data part
data = df.values

# Create the new file
with open('reformatted_output.csv', 'w') as f:
    # Write the ENSG IDs line
    f.write(ensg_line + '\n')
    
    # Write each row of data
    for row in data:
        # Convert each number to string and join with commas
        row_str = ','.join(str(x) for x in row)
        f.write(row_str + '\n')