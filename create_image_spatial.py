import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from collections import defaultdict
import os

def parse_gmt(gmt_file, gene_list):
    """
    Parse a GMT file to extract pathways and map genes to pathways.

    Parameters:
        gmt_file (str): Path to the GMT file.
        gene_list (list): List of gene symbols in the dataset (from expression file).

    Returns:
        dict: A dictionary where keys are pathway names and values are lists of gene indices.
    """
    pathway_map = defaultdict(list)
    with open(gmt_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            pathway_name = parts[0]
            genes = parts[2:]  # Skip the URL or description

            # Map to indices in the expression data
            gene_indices = [gene_list.index(gene) for gene in genes if gene in gene_list]

            if gene_indices:  # Only include pathways with mapped genes
                pathway_map[pathway_name] = gene_indices

    return pathway_map


def create_expression_images_with_cancer_genes(input_file, output_dir, gmt_file, target_size=(29, 29)):
    """
    Create images from cancer-related gene expression data, grouping genes by pathways.
    
    Parameters:
        input_file (str): Path to the gene expression CSV file (841 genes only).
        output_dir (str): Directory to save the images.
        gmt_file (str): Path to the GMT file with pathway definitions.
        target_size (tuple): Target size of the output images (rows, cols).
    """
    # Read the gene expression data
    data = pd.read_csv(input_file, index_col=0)

    # Get gene list and patient IDs
    gene_list = data.index.tolist()
    patient_ids = data.columns

    # Parse the GMT file and filter to only include genes in gene_list
    pathway_map = parse_gmt(gmt_file, gene_list)

    # Flatten the spatial map: Only include genes from the dataset grouped by pathways
    spatial_map = []
    for pathway, genes in pathway_map.items():
        spatial_map.extend(genes)

    # Check if the spatial map fits into the target size
    total_pixels = target_size[0] * target_size[1]
    if len(spatial_map) > total_pixels:
        spatial_map = spatial_map[:total_pixels]  # Truncate to fit the image size
    elif len(spatial_map) < total_pixels:
        spatial_map += [-1] * (total_pixels - len(spatial_map))  # Pad with -1

    # Find the global min and max values across all samples
    global_min = data.values.min()
    global_max = data.values.max()
    print(f"Global value range: [{global_min:.2f}, {global_max:.2f}]")

    # Get the values as a numpy array
    values = data.values

    # Calculate pixel values using the formula:
    def convert_to_pixel_values(x):
        pixel_values = np.round((x - global_min) * 255 / (global_max - global_min))
        return np.clip(pixel_values, 0, 255).astype(np.uint8)

    pixel_values = convert_to_pixel_values(values)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Create images for each sample
    for i, patient_id in enumerate(patient_ids):
        # Initialize the image array with zeros
        img_array = np.zeros(target_size, dtype=np.uint8)

        # Map the pixel values to the spatial arrangement
        for idx, gene_idx in enumerate(spatial_map):
            if gene_idx != -1:  # If the position is mapped to a gene
                row = idx // target_size[1]
                col = idx % target_size[1]
                img_array[row, col] = pixel_values[gene_idx, i]

        # Create and save image
        img = Image.fromarray(img_array)
        sample_output = os.path.join(output_dir, f'{patient_id}.png')
        img.save(sample_output)

        # Display first few images
        if i < 5:  # Show first 5 images
            plt.figure(figsize=(5, 5))
            plt.imshow(img_array, cmap='gray')
            plt.title(f'Patient {patient_id}\nPixel Range: [{np.min(img_array)}, {np.max(img_array)}]')
            plt.axis('off')
            plt.savefig(os.path.join(output_dir, f'preview_{patient_id}.png'))
            plt.close()

    print(f"Created images for {pixel_values.shape[1]} samples")
    print(f"Image shape: {target_size}")
    print(f"Images saved to: {output_dir}")
    print(f"Used global range [{global_min:.2f}, {global_max:.2f}] for scaling")


def parse_gmt(gmt_file, gene_list):
    """
    Parse a GMT file and return only pathways with genes in the gene_list.
    
    Parameters:
        gmt_file (str): Path to the GMT file.
        gene_list (list): List of gene symbols in the dataset.

    Returns:
        dict: A dictionary where keys are pathway names and values are lists of gene indices.
    """
    pathway_map = defaultdict(list)
    with open(gmt_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            pathway_name = parts[0]
            genes = parts[2:]  # Skip the URL or description
            # Only include genes that are in the provided gene_list
            gene_indices = [gene_list.index(gene) for gene in genes if gene in gene_list]
            if gene_indices:
                pathway_map[pathway_name] = gene_indices
    return pathway_map




if __name__ == "__main__":
    input_file = '/Users/stanleychen/git/Melanoma/final_data/RNA_cancer_related_genes.csv'  # Replace with your gene expression CSV file
    gmt_file = '/Users/stanleychen/git/Melanoma/clinical_data/data/gene_pathways/c2.cgp.v2024.1.Hs.symbols.gmt'  # Replace with your GMT file
    output_dir = 'image_data_spatial'
    target_size = (29, 29)
    max_pathways = 10  # Use only the top 10 pathways

    create_expression_images_with_cancer_genes(input_file, output_dir, gmt_file, target_size)
