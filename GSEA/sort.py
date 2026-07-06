import pandas as pd

# Path to your differential expression file
file_path = "/Users/stanleychen/git/Melanoma/GSEA/results/differential_expression.csv"

# Load the data
df = pd.read_csv(file_path, index_col=0)

# Ensure columns are numeric
df["padj"] = pd.to_numeric(df["padj"])

# Sort by padj values from lowest to highest
sorted_df = df.sort_values(by="padj")

# Display the top 20 most significant genes
print(sorted_df[["log2FoldChange", "pvalue", "padj"]].head(20))

# Optionally save the sorted data to a new file
sorted_df.to_csv("/Users/stanleychen/git/Melanoma/GSEA/results/differential_expression_sorted.csv")