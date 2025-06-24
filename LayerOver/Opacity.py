# -*- coding: utf-8 -*-
"""
Created:   2024-10-08 
Modified:  2024-10-08
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Wrapper class for opacity generating and analysis utilities.

"""
version = '0.1.0'
last_modified_date = '2024-10-08'

  #System and built-ins
import os
import sys
import json
import math
import csv
from tkinter import Tk, filedialog
from random import choice
import random

  #Visualizaiton
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

  #Data Handling
import numpy as np
import pandas as pd

  #Scientifiic algorithm packages
import scipy.interpolate
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

  #Utilities
from tqdm.auto import tqdm

  #Import other LayerUp modules
import Gcode



class Opacity:
   
    blank_spec_dict = {
        'n_voxel_points': 5,      #Number of points for that get evaluated for opacity in a part
        'n_pixels': 100,          #Number of pixels in 'Camera'; Number of pixels in final opacity measure
        'camera_radius': 1,       #Radius of points around voxel center to consider
        'distance_cutoff_threshold': 1,           #distance between pixel vector and layer point (in mm typically) 
        'max_opacity': 0.1,   # max decrement multiplier of light crossing through thickest part of strand
        }

    
    general_opacity_func = lambda distance: max_opacity * (math.sqrt(strand_radius**2 - distance**2) / strand_radius)
    
    
    def generate_opacity_from_part(Gcode_object = '', object_filenames = '', spec_dict = {}, num_points = 1):
    
        #Populate default values
        if len(spec_dict)<1:
            spec_dict = Opacity.blank_spec_dict
        
        for d_key in list(spec_dict.keys()):
            if 'n_voxel_points' in d_key:
                n_voxel_points = spec_dict[d_key]
            if 'n_pixels' in d_key:
                n_pixels = spec_dict[d_key]
            if 'camera_radius' in d_key:
                camera_radius = spec_dict[d_key]
            if 'distance_cutoff_threshold' in d_key:
                threshold = spec_dict[d_key]
            if 'max_opacity' in d_key:
                max_opacity = spec_dict[d_key]
                
        if 'opacity_func' not in list(spec_dict.keys()):
            opacity_func = Opacity.general_opacity_func
        
        if Gcode_object == '':
            if object_filenames =='':
                part = Gcode()
                part.get_gcode()
            else:
                part = Gcode()
                part.get_gcode(files=object_filenames)
        else:
            part = Gcode_object
        
        
          #get names of layers
        layer_list = list(part.layers.keys())
        
          #auto-generate part name guess and assign directory to save to
        part_dir = os.path.dirname(part.config_path)
        part_name_guess = os.path.basename(part_dir)
        
          #check if there's a single strand diameter
        if part.strand_diameter != 0:
            strand_radius = part.strand_diameter/2*10
            homogenous_strand_size = True
        else:
            homogenous_strand_size = False
            radius_list = [size/2*10 for size in part.layer_strand_diameters]

    ###############################################################################
    ####   Development functions   ################################################
    ###############################################################################

