"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   20YY-MM-DD
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: 

"""

from PIL import Image
import numpy as np
import csv
import os
from tkinter import filedialog, Tk

def create_grayscale_histogram(image_path, output_csv='histogram.csv'):
    """
    Read an image, create a grayscale histogram, and save to CSV.
    
    Args:
        image_path: Path to the input image
        output_csv: Path for the output CSV file (default: 'histogram.csv')
    """
    try:
        # Read the image using Pillow
        print(f"Reading image: {image_path}")
        img = Image.open(image_path)
        
        # Convert to grayscale if it's not already
        if img.mode != 'L':
            print("Converting image to grayscale...")
            img = img.convert('L')
        
        # Convert to numpy array
        img_array = np.array(img)
        
        # Create histogram using numpy
        # bins=256 for values 0-255, range=(0, 256) covers all grayscale values
        histogram, bin_edges = np.histogram(img_array.flatten(), bins=256, range=(0, 256))
        
        # Create dictionary with pixel value as key and count as value
        histogram_dict = {
            int(i): int(count) for i, count in enumerate(histogram)
        }
        
        # Save to CSV
        print(f"Saving histogram to: {output_csv}")
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write header
            writer.writerow(['Pixel_Value', 'Count'])
            # Write data
            for pixel_value, count in histogram_dict.items():
                writer.writerow([pixel_value, count])
        
        print(f"✓ Successfully created histogram with {len(histogram_dict)} bins")
        print(f"✓ Image dimensions: {img_array.shape}")
        print(f"✓ Total pixels: {img_array.size}")
        
        return histogram_dict
        
    except FileNotFoundError:
        print(f"Error: Image file '{image_path}' not found!")
        return None
    except Exception as e:
        print(f"Error: {str(e)}")
        return None

# Replace with your image path
root = Tk()
image_path = filedialog.askopenfilename(title="Select an image file")
root.destroy()
image_dir = os.path.dirname(image_path)
output_csv_filepath = os.path.join(image_dir, 'histogram.csv')

# Create histogram and save to CSV
histogram = create_grayscale_histogram(image_path, output_csv_filepath)
    
# Optional: Print some statistics
if histogram:
    print("\n--- Statistics ---")
    print(f"Min pixel value with count > 0: {min(k for k, v in histogram.items() if v > 0)}")
    print(f"Max pixel value with count > 0: {max(k for k, v in histogram.items() if v > 0)}")
    print(f"Most common pixel value: {max(histogram, key=histogram.get)} "
            f"(count: {histogram[max(histogram, key=histogram.get)]})")