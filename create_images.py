import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def create_expression_images(input_file, output_dir, target_size=(29, 29)):
    # Read the gene expression data
    data = pd.read_csv(input_file, index_col=0)
    
    # Get patient IDs from column names
    patient_ids = data.columns
    
    # Find the global min and max values across all samples
    global_min = data.values.min()
    global_max = data.values.max()
    print(f"Global value range: [{global_min:.2f}, {global_max:.2f}]")
    
    # Get the values as a numpy array
    values = data.values
    
    # Calculate pixel values using the formula:
    # PixelValue = Round((CellValue - MinValue) * 255 / (MaxValue - MinValue))
    def convert_to_pixel_values(x):
        pixel_values = np.round((x - global_min) * 255 / (global_max - global_min))
        return np.clip(pixel_values, 0, 255).astype(np.uint8)
    
    pixel_values = convert_to_pixel_values(values)
    
    # Calculate required padding
    total_pixels = target_size[0] * target_size[1]  # 29 * 29 = 841
    current_features = values.shape[0]  # 839
    padding_needed = total_pixels - current_features
    
    # Pad with zeros if needed
    if padding_needed > 0:
        padding = np.zeros((padding_needed, pixel_values.shape[1]), dtype=np.uint8)
        pixel_values = np.vstack([pixel_values, padding])
    
    # Create output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Create images for each sample
    for i, patient_id in enumerate(patient_ids):
        # Reshape to 29x29
        img_array = pixel_values[:, i].reshape(target_size)
        
        # Create and save image
        img = Image.fromarray(img_array)
        sample_output = os.path.join(output_dir, f'{patient_id}.png')
        img.save(sample_output)
        
        # Display first few images
        if i < 5:  # Show first 5 images
            plt.figure(figsize=(5, 5))
            plt.imshow(img_array, cmap='gray')
            plt.title(f'Patient {patient_id}\nValue Range: [{np.min(values[:,i]):.2f}, {np.max(values[:,i]):.2f}]\nPixel Range: [{np.min(img_array)}, {np.max(img_array)}]')
            plt.axis('off')
            plt.savefig(os.path.join(output_dir, f'preview_{patient_id}.png'))
            plt.close()
    
    print(f"Created images for {pixel_values.shape[1]} samples")
    print(f"Image shape: {target_size}")
    print(f"Original features: {current_features}")
    print(f"Padding added: {padding_needed}")
    print(f"Images saved to: {output_dir}")
    print(f"Used global range [{global_min:.2f}, {global_max:.2f}] for scaling")

if __name__ == "__main__":
    input_file = '/Users/stanleychen/git/Melanoma/final_data/RNA_cancer_related_genes.csv'
    output_dir = 'image_data'
    create_expression_images(input_file, output_dir)