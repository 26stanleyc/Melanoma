import re

def convert_to_tab_delimited(input_file, output_file):
    with open(input_file, 'r') as f:
        content = f.read()
        
    # Replace multiple spaces with a single tab
    formatted = re.sub(r'\s{2,}', '\t', content)
    
    with open(output_file, 'w') as f:
        f.write(formatted)

# Example usage
input_file = "/Users/stanleychen/git/Melanoma/R_scripts/GSEA_with_multiple_controls/filtered_ranked_list.tsv"
output_file = "/Users/stanleychen/git/Melanoma/R_scripts/GSEA_with_multiple_controls/filtered_ranked_list_v2.tsv"
convert_to_tab_delimited(input_file, output_file)