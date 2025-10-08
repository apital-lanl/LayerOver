# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 16:35:38 2025

@author: 359794
"""

import Opacity
from tkinter import Tk, filedialog

#Load filenames for each part
filenames_lists = []
continue_check = True
while continue_check:
    # Get the names by user selection
    root = Tk()
    these_filenames = filedialog.askopenfilenames(title= 'Select images to stitch', filetypes = [('LayerOver Files', '*.json *.pgm')])
    root.destroy()
    
    # Sort to get images in left-to-right, top-to-bottom order by index (hopefully)
    filenames_lists.append(these_filenames)
    
    # Ask if user want to add more
    user_report = input("Stitch another set of images (y/n or any key+ENTER to quit)").lower()
    if 'y' in user_report:
        continue_check = True
    else:
        continue_check = False


voxel_dict = Opacity.opacity_from_gcode(filenames_lists = filenames_lists, n_voxel_points = 5, 
                                n_pixels = 150, camera_radius = 3, 
                                threshold = 1, max_opacity= 0.1, 
                                opacity_func = 'fixed-radial', 
                                size_multiplier = 1,
                                strand_thickness_func = 'circular')