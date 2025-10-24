# -*- coding: utf-8 -*-
"""
? 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2024-10-08 
Modified:  2024-10-08
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Wrapper class for opacity generating and analysis utilities.

"""
version = '0.1.0'
last_modified_date = '2024-10-08'

  #System and built-ins
from operator import itemgetter
import os
import sys
import json
import math
import csv
from tkinter import Tk, filedialog
from random import choice
import random
import traceback

  #Visualizaiton
import matplotlib.pyplot as plt

  #Data Handling
import numpy as np
import pandas as pd

  #Utilities
from tqdm.auto import tqdm

  #Import other LayerUp modules
from Part import Gcode
import Points
import Voxel


blank_spec_dict = {
    'n_voxel_points': 5,            #Number of points for that get evaluated for opacity in a part
    'n_pixels': 100,                #Number of pixels in 'Camera'; Number of pixels in final opacity measure
    'camera_radius': 1,             #Radius of points around voxel center to consider
    'distance_cutoff_threshold': 1, #distance between pixel vector and layer point (in mm typically) 
    'max_opacity': 0.1,             #decrease/250um; max decrement multiplier of light crossing through thickest part of strand
    }
general_opacity_func = lambda distance: max_opacity * (math.sqrt(strand_radius**2 - distance**2) / strand_radius)
default_full_strand_opacity = 0.1           #decrease/250um; max decrement multiplier of light crossing through thickest part of strand
default_opacity_function = 'fixed-radial'   #
default_squish_factor= 0.3                  #compaction ratio vs. ideal stack (0.3 = 30% compaction, i.e. 1" part ends up being 0.7")
default_base_compaction_ratio= 0.2          #ratio of total 'squish_factor' that occurs at the baseplate; i.e. '0.3' means 30% of normal compaction occurs at baseplate

def opacity_from_gcode(filenames_lists = [], n_voxel_points = 5, n_pixels = 100, camera_radius = 3,
                    max_opacity= default_full_strand_opacity, 
                    opacity_func = default_opacity_function,
                    threshold = 1, size_multiplier = 1, 
                    strand_thickness_func = 'circular'):
    
    ''' v1.0.1  created: 2024-10-06  modified:  2025-07-28
    
    Description: Take layer gcode files and a configuration file (LayerUp style) and get opacity measures.
        Can handle multiple parts in succession.
    
    INPUT:
        filenames_lists     List of lists; [ ['config', 'layer1', ...],]    
          #voxel/camera parameters
        n_voxel_points      Number of points for that get evaluated for opacity in a part
        n_pixels            Number of pixels in 'Camera'; Number of pixels in final opacity measure
        camera_radius       Radius of points around voxel center to consider
          
          #Opacity measure parameters
        threshold           Distance between pixel vector and layer point (in mm typically) 
        max_opacity         Max decrement multiplier of light crossing through thickest part of strand
        opacity_func        Function to use for calculating opacity from distance-to-strand-center normal measure
            'fixed-radial'- assumes a constant, round shape 
    '''

    squish_factor= default_squish_factor                    #compaction ratio vs. ideal stack (0.3 = 30% compaction, i.e. 1" part ends up being 0.7")
    base_compaction_ratio= default_base_compaction_ratio    #ratio of total 'squish_factor' that occurs at the baseplate; i.e. '0.3' means 30% of normal compaction occurs at baseplate
    
    if len(filenames_lists) < 1:
        continue_check = True
        filenames_lists = []
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
    
    #TODO: add options in the future
      #if 'strand_thickness_func' is already a callable object, assume it's a function and just use it
    if callable(strand_thickness_func):
        pass
    #otherwise, start parsing what the input wants the functional form to look like
      #standard circular strand assumptions
    elif strand_thickness_func == 'circular':
        strand_thickness_func = lambda distance: 2* math.sqrt(strand_radius**2 - distance**2)

    #TODO: add options in the future
    if opacity_func == 'fixed-radial':
        opacity_func = lambda distance: max_opacity * strand_thickness_func(distance) / (2*strand_radius) 
     
    #Default to 'fixed-radial'
    else:
        opacity_func = lambda distance: max_opacity * strand_thickness_func(distance) / (2*strand_radius)
    
    #%% (3) Do everything
    for file_list in filenames_lists:
        try:
            part = Gcode()
            part.get_gcode(files = file_list)
            
              #get names of layers
            layer_list = list(part.layers.keys())
            
              #auto-generate part name guess and assign directory to save to
            part_dir = os.path.dirname(part.config_path)
            part_name_guess = os.path.basename(part_dir)
            
              #check if there's a single strand diameter
            if part.strand_diameter != 0:
                strand_radius = part.strand_diameter/2*size_multiplier
                homogenous_strand_size = True
                radius_list = [size/2*size_multiplier for size in part.layer_strand_diameters]
            else:
                homogenous_strand_size = False
                radius_list = [size/2*size_multiplier for size in part.layer_strand_diameters]
                
            #Generate center point, normals, etc. from fibonacci points and near neighbors
            point_dicts = Voxel.get_random_point_summaries(part, n_retrieval_points = n_voxel_points, \
                                               excluded_point_indcs = [], show_iterations = True)

            # Make voxels and visualize
            voxel_opacity_images = []
            for idx, key in enumerate(list(point_dicts.keys())):
                
                #generate voxel name
                voxel_name = str(part_name_guess + "_voxel-" + str(idx))
                voxel_point_name = str(part_name_guess + "_VoxelPoints" + str(idx))
                save_name = os.path.join(part_dir, voxel_point_name)
                
                #pull voxel features
                center = point_dicts[idx]['approximate_center']
                normal_vector = point_dicts[idx]['unit_normal']
                radius = camera_radius      #TODO: leave option to make this layer-defined
                #TODO: Radius is currently hard-coded by user or default input; could be made to be more flexible here
                
                #get voxel bounds and add to dictionary 
                voxel_dict = Voxel.get_voxel(center, normal_vector, radius= radius, voxel_depth= 2*radius, n_pixels=n_pixels)
                voxel_dict = Voxel.get_layer_points(part, voxel_dict, boundary_buffer_multiplier = 0.2)

                #update voxel_dict
                voxel_dict.update({'layer_names': layer_list})
                voxel_dict.update({'voxel_name': voxel_name})
                voxel_dict.update({'save_filepath': save_name})
                
                #Visualize voxel
                Voxel.voxel_points_visualize(voxel_dict, 
                           gif = False,
                           gif_frames = 4,
                           gif_integrals = 50,
                           azim_total_degrees = 360,
                           azim_step = 0,
                           save_plot = True)
                    
                #Flatten array of coordinates to make it homogenous
                #   (otherwise it's a list of lists, and each sub-list is a different size; this makes numpy sad)
                voxel_dict.update( {'flat_layer_coordinates':{} } )
                for layer_key in list(part.layers.keys()):
                    #flatten coordinate array; kept previously as inhomogenous list to make it easier to parse consecutive strand coordinates
                       #NOTE: each 'this_list' is a strand of the print
                    flat_coordinates = []
                    for this_list in voxel_dict['layer_coordinates'][layer_key]['coordinates']:
                        for point in this_list:
                            flat_coordinates.append(point)
                
                    #re-wrap these lists as numpy arrays to make the following stuff faster and the syntax clearer (mostly the faster part)
                    flat_coordinates = np.array(flat_coordinates)
                    flat_indices = np.array(voxel_dict['layer_coordinates'][layer_key]['indices'])
                
                    voxel_dict['flat_layer_coordinates'].update ( {layer_key: {
                                            'coordinates':flat_coordinates, \
                                            'indices':flat_indices} })
                
                #Get the actual distances for each pixel and calculate opacity
                #NOTE: this is v1 and isn't vectorized, making it incredibly slow
                pixel_min_distances = []
                layer_opacity_images = []
                for l_idx, layer_key in enumerate(list(part.layers.keys())):
                
                    flat_coordinates = voxel_dict['flat_layer_coordinates'][layer_key]['coordinates']
                    flat_indices = voxel_dict['flat_layer_coordinates'][layer_key]['indices']
                    layer_opacity_array = np.ones( (voxel_dict['top_pixel_coordinate_array'].shape[0]))
                
                    for pixel_idx in tqdm(range(voxel_dict['top_pixel_coordinate_array'].shape[0])):
                        top_pixel = voxel_dict['top_pixel_coordinate_array'][pixel_idx]
                        bot_pixel = voxel_dict['bottom_pixel_coordinate_array'][pixel_idx]
                        
                        #downselect closest points
                        dist_func = lambda pbot, ptop, p3: np.linalg.norm(np.cross(pbot-ptop, ptop-p3))/np.linalg.norm(pbot-ptop)
                        coord_distances = np.array([dist_func(bot_pixel, top_pixel, layer_point) for layer_point in flat_coordinates])
                        close_mask = coord_distances <= threshold   #get boolean array of 'close' layer points
                        close_indices = flat_indices[close_mask]
                
                        #for each "close point" index, get previous and next point indices (if they exist) and add them to a running list
                        close_segment_points = []
                        skip_last_index = False  #initialize
                        skip_next_index = False  #initialize
                        for this_point_idx in close_indices:
                            this_point = part.layers[layer_key]['coordinates'][this_point_idx]
                            
                            next_point_idx = this_point_idx+1
                            last_point_idx = this_point_idx-1
                            try:
                                last_point = part.layers[layer_key]['coordinates'][last_point_idx]
                            except:
                                skip_last_index = True
                            
                            try:
                                next_point = part.layers[layer_key]['coordinates'][next_point_idx]
                            except:
                                skip_next_index = True
                        
                            #make sure each 
                            if (not skip_next_index) and (not skip_last_index):
                                close_segment_points.append( [last_point, this_point] )
                                close_segment_points.append( [this_point, next_point] )
                                
                            elif skip_next_index:
                                close_segment_points.append( [last_point, this_point] )
                            
                            elif skip_last_index:
                                close_segment_points.append( [this_point, next_point] )
                        close_segment_points = np.array(close_segment_points)
                        
                        #get minimum segment-to-segment perpendicular distance
                        distances = np.array([Points.line_line_distance(layerp1, layerp2, bot_pixel, top_pixel)[2]
                                     for layerp1, layerp2 in close_segment_points])
                        try:
                            min_distance = np.min(distances)
                        #if 'distances' is empty, just set to minimum opacity distance
                        #TODO: this is an inelegant solution; need better scoping here
                        except ValueError:
                            min_distance = strand_radius + (strand_radius*0.1)
                        #check if strand radius is homogenous; if so, use that value, otherwise get layer's individual radius
                        #TODO: add difference flags if the layers have different strand diameters
                        if homogenous_strand_size:
                            strand_radius = radius_list[l_idx]
                            max_opacity = max_opacity * (strand_radius/250) #adjust max opacity to account for 250um nozzle nominal assumption
                        else:
                            strand_radius = radius_list[l_idx]
                            max_opacity = max_opacity * (strand_radius/250) #adjust max opacity to account for 250um nozzle nominal assumption
                        
                          #get an opacity if there's actually a strand in this pixel
                        if min_distance <= strand_radius:
                            this_opacity = opacity_func(min_distance)
                        else:
                            this_opacity = 1
                        
                          #update opacity image
                        # opacity_array[pixel_idx] = opacity_array[pixel_idx] * this_opacity
                        layer_opacity_array[pixel_idx] = this_opacity

                    layer_opacity_image = np.reshape(layer_opacity_array, (n_pixels, n_pixels))
                    layer_opacity_images.append(layer_opacity_image)
                    plt.imshow(layer_opacity_image)
                    plt.title(f"{voxel_name}-{layer_key}")
                    save_name = os.path.join(part_dir, str(f"{voxel_name}-{layer_key}"+ '.png'))
                    plt.imsave(save_name, layer_opacity_image, dpi=600)
                    plt.show()

                    partial_layer_name = f"_{voxel_name}-{layer_key}.npy"
                    save_name = os.path.join(part_dir, partial_layer_name)
                    np.save(save_name, layer_opacity_image)
                
                #Combine each layer opacity into a single voxel's opacity
                #  Adjusts for the fact that strands coallesce to some extent
                for layer_idx, layer_image in enumerate(layer_opacity_images):
                    #Make adjustments for f'd up base layer
                    if layer_idx == 0:
                        max_squished_opacity = max_opacity * (1-(squish_factor*base_compaction_ratio))  #assume baselayer squishes all over at max-strand-thickness, by about 'base_compaction_ratio'
                        layer_image[layer_image<max_squished_opacity]= max_squished_opacity
                        opacity_array = layer_image
                    #Make adjustments specifically for second layer (interacts w/ f'd up base layer)
                    elif layer_idx == 1:
                        sum_image = layer_opacity_images[layer_idx-1] * layer_image
                        max_squished_opacity = max_opacity**2 * (1-squish_factor)  #Each overlapping region will be (1-squish_factor) less opaque than expected 
                        max_layer_opacity = max_opacity * (1-squish_factor)
                        layer_image[sum_image<max_squished_opacity] = max_layer_opacity
                        opacity_array = opacity_array * layer_image

                    else:
                        sum_image = layer_opacity_images[layer_idx-1] * layer_image
                        max_squished_opacity = max_opacity**2 * (1-squish_factor)  #Each overlapping region will be (1-squish_factor) less opaque than expected 
                        max_layer_opacity = max_opacity * (1-squish_factor)
                        layer_image[sum_image<max_squished_opacity] = max_layer_opacity
                        opacity_array = opacity_array * layer_image
                
                #Clean up the opacity image
                opacity_image = np.reshape(opacity_array, (n_pixels, n_pixels))
                voxel_opacity_images.append(opacity_image)
                #Save the numpy array
                partial_name = f"_{voxel_name}_CombinedOpacity.npy"
                save_name = os.path.join(part_dir, partial_name)
                np.save(save_name, opacity_image)
                #Show and save the opacity image
                plt.imshow(opacity_image)
                plt.title(f"{voxel_name}-CombinedOpacity")
                partial_name =  f"_{voxel_name}_CombinedOpacity.png"
                save_name = os.path.join(part_dir, partial_name)
                plt.imsave(save_name, opacity_image, dpi=600)
                plt.show()
            
            opacity_arrays = np.array(voxel_opacity_images)
            partial_name = str(part_name_guess + '_OpacityArrays.npy')
            save_name = os.path.join(part_dir, partial_name)
            np.save(save_name, opacity_arrays)
        except IndexError as IE:
            print()
            print(IE)
            print()
            print(traceback.format_exc())   
            print()

    voxel_dict.update({'point_dicts':point_dicts})
    voxel_dict.update({'radius_list': radius_list})

    return voxel_dict


def opacity_from_voxel_dict(voxel_dict = None, file_list = None, n_voxel_points = 1, n_pixels = 100, camera_radius = 3,
                    threshold = 1, max_opacity= 0.1, opacity_func = 'fixed-radial', size_multiplier = 1,
                    strand_thickness_func = 'circular'):
    
    ''' v0.1.0  created: 2025-07-28  modified:  2025-07-28

    NOTE: Not a full alpha; needs substantial work
    
    Description: Get opacity straight from a 'voxel_dict' that includes flattened voxel coordinate lists.
        If no flattened coordinates exist, try to make them. If that fails, make a new voxel_dict from user-defined
        files and use that voxel dict instead
    
    INPUT:
        file_list           List containing at least one layer (gcode) file and a config.json file
          #voxel/camera parameters
        n_voxel_points      Number of points for that get evaluated for opacity in a part
        n_pixels            Number of pixels in 'Camera'; Number of pixels in final opacity measure
        camera_radius       Radius of points around voxel center to consider
          #Opacity measure parameters
        threshold           Distance between pixel vector and layer point (in mm typically) 
        max_opacity         Max decrement multiplier of light crossing through thickest part of strand
        opacity_func        Function to use for calculating opacity from distance-to-strand-center normal measure
            'fixed-radial'- assumes a constant, round shape

    TODO:
        - fix first pass alpha results
    '''

    squish_factor = 0.3    #compaction ratio vs. ideal stack (0.3 = 30% compaction, i.e. 1" part ends up being 0.7")

    if voxel_dict == None:
        if file_list == None:
            root = Tk()
            file_list = filedialog.askopenfilenames(title= 'Select files for a single Part', filetypes=[("LayerUp Files", "*.json *.pgm")])
            root.destroy()

        voxel_dict = Voxel.get_voxels_from_filelist(file_list, n_voxel_points = n_voxel_points, 
                            n_pixels = n_pixels, camera_radius = camera_radius,
                            size_multiplier =size_multiplier)
    
    #TODO: add options in the future
      #if 'strand_thickness_func' is already a callable object, assume it's a function and just use it
    if callable(strand_thickness_func):
        pass
    #otherwise, start parsing what the input wants the functional form to look like
      #standard circular strand assumptions
    elif strand_thickness_func == 'circular':
        strand_thickness_func = lambda distance: 2* math.sqrt(strand_radius**2 - distance**2)

    #TODO: add options in the future
    if opacity_func == 'fixed-radial':
        opacity_func = lambda distance: max_opacity * strand_thickness_func(distance) / (2*strand_radius) 
     
    #Default to 'fixed-radial'
    else:
        opacity_func = lambda distance: max_opacity * strand_thickness_func(distance) / (2*strand_radius)
    
    point_dicts = voxel_dict['point_dicts']
    radius_list = voxel_dict['radius_list']
    save_name = voxel_dict['save_filepath']
    voxel_name = voxel_dict['voxel_name']
    voxel_point_name = os.path.basename(save_name)
    part_dir = os.path.dirname(save_name)
    part_name_guess = os.path.basename(part_dir)
    layer_keys = list(voxel_dict['flat_layer_coordinates'])
    radius_list = voxel_dict['radius_list']
        
    #%% (3) Do everything
    try:
        # Make voxels and visualize
        voxel_opacity_images = []
        for idx, key in enumerate(list(point_dicts.keys())):
                
            #Visualize voxel
            Voxel.voxel_points_visualize(voxel_dict, 
                        gif = False,
                        gif_frames = 4,
                        gif_integrals = 50,
                        azim_total_degrees = 360,
                        azim_step = 0,
                        save_plot = True)
                
            #Get the actual distances for each pixel and calculate opacity
            #NOTE: this is v1 and isn't vectorized, making it incredibly slow
            pixel_min_distances = []
            opacity_array = np.ones( (voxel_dict['top_pixel_coordinate_array'].shape[0]))
            for l_idx, layer_key in enumerate(layer_keys):
                
                flat_coordinates = voxel_dict['flat_layer_coordinates'][layer_key]['coordinates']
                flat_indices = voxel_dict['flat_layer_coordinates'][layer_key]['indices']
                layer_opacity_array = np.ones( (voxel_dict['top_pixel_coordinate_array'].shape[0]))
                
                for pixel_idx in tqdm(range(voxel_dict['top_pixel_coordinate_array'].shape[0])):
                    top_pixel = voxel_dict['top_pixel_coordinate_array'][pixel_idx]
                    bot_pixel = voxel_dict['bottom_pixel_coordinate_array'][pixel_idx]
                        
                    #downselect closest points
                    dist_func = lambda pbot, ptop, p3: np.linalg.norm(np.cross(pbot-ptop, ptop-p3))/np.linalg.norm(pbot-ptop)
                    coord_distances = np.array([dist_func(bot_pixel, top_pixel, layer_point) for layer_point in flat_coordinates])
                    close_mask = coord_distances <= threshold   #get boolean array of 'close' layer points
                    close_indices = flat_indices[close_mask]
                
                    #for each "close point" index, get previous and next point indices (if they exist) and add them to a running list
                    close_segment_points = []
                    skip_last_index = False  #initialize
                    skip_next_index = False  #initialize
                    for this_point_idx in close_indices:
                        this_point = flat_coordinates[this_point_idx]
                            
                        next_point_idx = this_point_idx+1
                        last_point_idx = this_point_idx-1
                        try:
                            last_point = flat_coordinates[last_point_idx]
                        except:
                            skip_last_index = True
                            
                        try:
                            next_point = flat_coordinates[next_point_idx]
                        except:
                            skip_next_index = True
                        
                        #make sure each 
                        if (not skip_next_index) and (not skip_last_index):
                            close_segment_points.append( [last_point, this_point] )
                            close_segment_points.append( [this_point, next_point] )
                                
                        elif skip_next_index:
                            close_segment_points.append( [last_point, this_point] )
                            
                        elif skip_last_index:
                            close_segment_points.append( [this_point, next_point] )
                    close_segment_points = np.array(close_segment_points)
                        
                    #get minimum segment-to-segment perpendicular distance
                    distances = np.array([Points.line_line_distance(layerp1, layerp2, bot_pixel, top_pixel)[2]
                                    for layerp1, layerp2 in close_segment_points])
                    try:
                        min_distance = np.min(distances)
                    #if 'distances' is empty, just set to minimum opacity distance
                    #TODO: this is an inelegant solution; need better scoping here
                    except ValueError:
                        min_distance = strand_radius + (strand_radius*0.1)
                    
                    #Pull strand radius (hopefully) from 'voxel_dict' radius list
                    strand_radius = radius_list[l_idx]
                        
                        #get an opacity if there's actually a strand in this pixel
                    if min_distance <= strand_radius:
                        this_opacity = opacity_func(min_distance)
                    else:
                        this_opacity = 1
                        
                        #update opacity image
                    opacity_array[pixel_idx] = opacity_array[pixel_idx] * this_opacity
                    layer_opacity_array[pixel_idx] = this_opacity

                layer_opacity_image = np.reshape(layer_opacity_array, (n_pixels, n_pixels))
                plt.imshow(layer_opacity_image)
                plt.title(f"{voxel_name}-{layer_key}")
                save_name = os.path.join(part_dir, str(f"{voxel_name}-{layer_key}"+ '.png'))
                plt.imsave(save_name, layer_opacity_image, dpi=600)
                plt.show()

                partial_layer_name = str(part_name_guess + f"_{voxel_name}-{layer_key}.npy")
                save_name = os.path.join(part_dir, partial_layer_name)
                np.save(save_name, layer_opacity_image)
                    
            opacity_image = np.reshape(opacity_array, (n_pixels, n_pixels))
            voxel_opacity_images.append(opacity_image)
                
            plt.imshow(opacity_image)
            plt.title(f"{voxel_name}-CombinedOpacity")
            save_name = os.path.join(part_dir, str(voxel_name+ '.png'))
            plt.imsave(save_name, opacity_image, dpi=600)
            plt.show()
            
        opacity_arrays = np.array(voxel_opacity_images)
        partial_name = str(part_name_guess + '_OpacityArrays.npy')
        save_name = os.path.join(part_dir, partial_name)
        np.save(save_name, opacity_arrays)
    except IndexError as IE:
        print()
        print(IE)
        print()
        print(traceback.format_exc())   
        print()

    return voxel_dict


def adjust_opacity_for_squish(layer_arrays,
                              combined_opacity_array = None,
                              layer_index = None,
                              max_opacities = None,
                              full_strand_opacity = default_full_strand_opacity,
                              opacity_function = default_opacity_function,
                              squish_factor = default_squish_factor,
                              compaction_ratio = default_base_compaction_ratio):
    
    ''' v0.1.0  created: 2025-08-05  modified:  2025-08-06
    
    Description: Run through opacity arrays for each layer of a print and combine them based on 'squishy' factors. 
        I.e. strands coallesce with each other and with the underlying substrate, so modify the opacities to reflect that loss of material.
    
    INPUT:
        layer_arrays-      iterable of arrays for each layer; can handle lists or numpy arrays

    TODO:
        - lor
    '''

    #Make sure 'layer_arrays' is a numpy array and define its shape
    if type(layer_arrays) == list:
        layer_arrays = np.array(layer_arrays)
    layer_arrays_shape = layer_arrays.shape
    layer_shape = tuple(layer_arrays_shape[1::])

    if combined_opacity_array == None:
        opacity_array = np.ones(layer_shape)

    #If no layer defined, run the full array
    if layer_index == None:
        layer_index = layer_arrays_shape[0]-1  #Assume first index is the number of layers; otherwise wtf is this array, you know?
        indices_list = list(range(layer_index+1))
    #If select layers are defined, pretend first layer is a base layer and just run with it
    #TODO: have a think about whether this should/could do something useful
    elif (type(layer_index) == list):
        indices_list = layer_index
    else:
        indices_list = list(range(layer_index+1))
    
    #If no 'max_opacity' has been passed, just use the layer values to find a (hopefully) appropriate one
    if max_opacities == None:
        max_opacities = []
        for layer_index in indices_list:
            this_max = layer_arrays[layer_index, ::].max()
            max_opacities.append(this_max)
    #If a list of 'max_opacities' is already defined, try and use them.
    elif type(max_opacities) == list:
        if len(max_opacities)==len(indices_list):
            pass #great
        else:
            pass #TODO: handle this case, but it's an edge case and will error out right away
    #If the max opacitiy is a float or int, just use that for everything
    #TODO: add case that includes numpy numerical types
    elif (type(max_opacities)==int) or (type(max_opacities)==float):
        max_opacities = [max_opacities for idx in (range(len(indices_list)+1))]
    
    #Combine each layer opacity into a single voxel's opacity
    for layer_idx, max_opacity in zip(indices_list, max_opacities):

        layer_image = layer_arrays[layer_idx, ::]

        #Make adjustments for f'd up base layer
        if layer_idx == 0:
            max_squished_opacity = max_opacity * (1-(squish_factor * compaction_ratio))  #assume baselayer squishes all over at max-strand-thickness, by about 'base_compaction_ratio'
            layer_image[layer_image<max_squished_opacity]= max_squished_opacity
            opacity_array = layer_image
        #Make adjustments specifically for second layer (interacts w/ f'd up base layer)
        elif layer_idx == 1:
            sum_image = layer_arrays[layer_idx-1, ::] * layer_image
            max_squished_opacity = max_opacity**2 * (1-squish_factor)  #Each overlapping region will be (1-squish_factor) less opaque than expected 
            max_layer_opacity = max_opacity * (1-squish_factor)
            layer_image[sum_image<max_squished_opacity] = max_layer_opacity
            opacity_array = opacity_array * layer_image
        else:
            sum_image = layer_arrays[layer_idx-1] * layer_image
            max_squished_opacity = max_opacity**2 * (1-squish_factor)  #Each overlapping region will be (1-squish_factor) less opaque than expected 
            max_layer_opacity = max_opacity * (1-squish_factor)
            layer_image[sum_image<max_squished_opacity] = max_layer_opacity
            opacity_array = opacity_array * layer_image
    
    return opacity_array


def process_directory_numpy_arrays(directory=None,
                                   save_combined_json = True,
                                   save_combined_csv = True):
    
    ''' v0.1.0  created: 2025-08-05  modified:  2025-08-05
    
    Description: Walk a selected directory, identify all numpy arrays associated with toolpaths for voxel layers, 
        and create a new combined opacity and associated dictionary of summary stats.
    
    INPUT:
        directory-      if empty, open a dialog to select a directory to walk

    TODO:
        - lor
    '''

    #Initialize variables
    opacity_dict = {}
    numpy_array_filenames = []
    numpy_array_filepaths = []

    #Select directory if none provided or a non-path is provided
    if directory == None:
        root = Tk()
        directory = filedialog.askdirectory(title="Select a directory with numpy layer arrays")
        root.destroy()
    elif not os.isdir(directory):
        root = Tk()
        directory = filedialog.askdirectory(title="Select a directory with numpy layer arrays")
        root.destroy()

    #Walk the directory and grab all numpy arrays
    for root, dirs, files in os.walk(directory, topdown = True):
        for file in files:
            this_filepath = os.path.join(root, file)
            immediate_directory_basename = os.path.basename(root)
        
            #While walking, look for SEAM filetypes
            main_filetype = file.split('.')[-1].lower()
            if (main_filetype == "npy"):
                #Do the filename stuff
                this_name = file.replace('.npy', '').lower()
                this_filepath = os.path.join(root, file)
                #Make sure only the layer files are being added; otherwise things get annoying
                if ('voxel' in this_name) and ('layer' in this_name):
                    numpy_array_filenames.append(this_name)
                    numpy_array_filepaths.append(this_filepath)
   
    #Run through collected numpy arrays and partition results
    for this_name, this_filepath in zip(numpy_array_filenames, numpy_array_filepaths):
        name_split_list= this_name.split('_')
        #Assuming structure like "5H-01_5H-01_voxel-X_layerX"
        structure_name = name_split_list[0]
        voxel_name = 'unk'
        layer_name = 'unk'
        for part in name_split_list:
            if 'voxel' in part.lower():
                voxel_name = part
            if 'layer' in part.lower():
                layer_name = part
            #Covers case if underscores have been naughty
            if voxel_name == layer_name:
                split_name = voxel_name.split('-')
                voxel_idx = 10  #dummy variable to keep unboundlocalerror happy
                for idx, part in enumerate(split_name):
                    if 'voxel' in part:
                        voxel_name = part
                        voxel_idx = idx
                    if 'layer' in part:
                        layer_name = part
                    if idx == (voxel_idx +1):
                        voxel_name = f"{voxel_name}-{part}"


        #Try adding to dict, update dict if it doesn't work
        try:
            name_list = opacity_dict[structure_name]['numpy_filenames']
            name_list.append(this_name)
            opacity_dict[structure_name]['numpy_filenames'] = name_list
            filepath_list = opacity_dict[structure_name]['numpy_filepaths']
            filepath_list.append(this_filepath)
            opacity_dict[structure_name]['numpy_filepaths'] = filepath_list
            
            #Try adding to existing voxel, otherwise add voxel to structure dict
            try:
                name_list = opacity_dict[structure_name][voxel_name]['layer_names']
                name_list.append(layer_name)
                opacity_dict[structure_name][voxel_name]['layer_names'] = name_list
                filepath_list = opacity_dict[structure_name][voxel_name]['layer_filepaths']
                filepath_list.append(this_filepath)
                opacity_dict[structure_name][voxel_name]['layer_filepaths'] = filepath_list
            except KeyError:
                opacity_dict[structure_name].update({
                    voxel_name:{
                        'layer_names': [layer_name],
                        'layer_filepaths': [this_filepath],
                        'layer_stats':{},
                        'combined_opacity_stats':{}
                        }
                    }
                                                    )
        except KeyError:
            opacity_dict.update({
                structure_name: {
                    'numpy_filenames': [this_name],
                    'numpy_filepaths': [this_filepath],
                    voxel_name:{
                        'layer_names': [layer_name],
                        'layer_filepaths': [this_filepath],
                        'layer_stats':{},
                        'combined_opacity_stats':{}
                        }
                    }
                }
                                )

    #Run through partitioned results and generate combined opacity
    for structure_name in list(opacity_dict.keys()):
        
        key_list = [this_key for this_key in list(opacity_dict[structure_name].keys()) if ('numpy' not in this_key)]
        for voxel_name in key_list:
            voxel_dict = opacity_dict[structure_name][voxel_name]
            layer_names = voxel_dict['layer_names']
            layer_filepaths = voxel_dict['layer_filepaths']
            layer_indices = []
            #Convoluted solution, but accomodates existing layer gcode files and is mildly flexible
            #TODO: add support for other formats of layer indexing (i.e. "layer-XXX", "layer(X)", etc.)
            for name, filepath in zip(layer_names, layer_filepaths):
                last_numerical_idx = 0 
                for check_idx in range(-1,-4,-1):
                    if (name[check_idx].isnumeric()) and (last_numerical_idx > check_idx):
                        last_numerical_idx = check_idx
                layer_idx = ''
                for idx in range(last_numerical_idx, 0, 1):
                    layer_idx = str(layer_idx) + str(name[idx])
                layer_idx = int(layer_idx)
                this_tuple = (name, layer_idx, filepath)
                layer_indices.append(this_tuple)
    
            sorted_layers = sorted(layer_indices, key=itemgetter(1))
            layers_array = []
            for layer_tuple in sorted_layers:
                this_layer_name = layer_tuple[0]
                #Open the numpy array
                this_filepath = layer_tuple[2]
                this_layer_array = np.load(this_filepath)
                layers_array.append(this_layer_array)
                try:
                    #Run some stats on this layer's opacity alone and add to dictionary
                    this_summary_dict = get_opacity_array_statistics(this_layer_array)
                    voxel_dict['layer_stats'].update({
                        this_layer_name: this_summary_dict
                        })
                #Flagged when an array is zero-size
                #TODO: find out why array is zero-sized... shouldn't matter much for the following results
                except ValueError as VE:
                    print('_'*50)
                    print(f"Error on {structure_name}-{voxel_name}")
                    print(f"\t {VE}")
                    print()
                    print(f"\t {traceback.format_exc()}") 
                    print()
                
            #Add all the layers together, hopefully accounting for 'squish' factors
            new_layers_array = adjust_opacity_for_squish(layers_array)
            
            try:
                #Run some stats on this layer's opacity alone and add to dictionary
                this_summary_dict = get_opacity_array_statistics(new_layers_array)
                voxel_dict['combined_opacity_stats'].update(this_summary_dict)
                
                #Save the combined array
                part_dir = os.path.dirname(this_filepath)
                partial_name = f"{structure_name.upper()}_{voxel_name}_CombinedOpacity.npy"
                save_name = os.path.join(part_dir, partial_name)
                np.save(save_name, new_layers_array)
                
                #Save the 'voxel_dict' as a JSON file
                partial_name = f"{structure_name.upper()}_{voxel_name}_OpacityDictionary.json"
                save_name = os.path.join(part_dir, partial_name)
                with open(save_name, "w") as filename:
                    json.dump(voxel_dict, filename, indent=4)
            except ValueError as VE:
                print('_'*50)
                print(f"Error on {structure_name}-Full Combined Array")
                print(f"\t {VE}")
                print()
                print(f"\t {traceback.format_exc()}") 
                print()
                
    #Convert the dictionary to a pandas.DataFrame for outputting
    for structure_idx, structure_name in enumerate(list(opacity_dict.keys())):
        voxel_names = [key for key in list(opacity_dict[structure_name].keys()) if ('numpy' not in key.lower())]
        for voxel_idx, voxel_name in enumerate(voxel_names):
            #Pull individual voxel's combined opacity stats
            this_dict = opacity_dict[structure_name][voxel_name]['combined_opacity_stats']
            this_dict.update({'name': f"{structure_name}-{voxel_name}"})
            bad_columns = [column for column in list(this_dict.keys()) if ('histogram' in column.lower())]
            for bad_column in bad_columns:
                this_dict.pop(bad_column)
            
            #If it's the first iteration, create a DataFrame with appropriate columns
            if (voxel_idx == 0) and (structure_idx == 0):
                opacity_df = pd.DataFrame(this_dict, index = [0])
            else:
                new_df = pd.DataFrame(this_dict, index = [0])
                opacity_df = pd.concat([new_df, opacity_df], axis=0)
                
    if save_combined_json:
        #Save the full 'opacity_dict' as a JSON file
        partial_name = f"{structure_name.upper()}_AllVoxels_CombinedOpacityDictionary.json"
        save_name = os.path.join(directory, partial_name)
        with open(save_name, "w") as filename:
            json.dump(opacity_dict, filename, indent=4)
            
    if save_combined_csv:
        partial_name = f"{structure_name.upper()}_AllVoxels_CombinedOpacityMetrics.csv"
        save_name = os.path.join(directory, partial_name)
        opacity_df.reset_index(drop=True, inplace=True)
        with open(save_name, "w") as filename:
            opacity_df.to_csv(filename, index=False, lineterminator='\n')

    return opacity_df
        
        
def get_opacity_array_statistics(array,
                                 n_histogram_bins = 50,
                                 n_rounding_places= 4):
    ''' v0.1.0  created: 2025-08-06  modified:  2025-08-06
    
    Description: Standardized reporting function for an array of opacities.
    
    INPUT:
        array-  Array or numpy.array of opacities for a voxel or voxel-like region

    TODO:
        - lor
    '''
    
    #Process input and initialize variables
    summary_dict = {
        'max_opacity': 0,
        'min_opacity': 0,
        'stdev_opacity': 0,
        'average_opacity': 0,
        'opacity_sum': 0,
        'opaque_region_sum': 0,
        'opaque_fraction': 0,
        'transmittance': 0,
        'absorbance': 0,
        'histogram_bins': [],
        'histogram_cnts': [],
        }
    array = np.array(array)
    
    #NOTE: converting everything back to python native formats to allow for exporting dictionary as a JSON; otherwise it will shit itself
    #Max, min, std dev., avg, sum
    summary_dict['max_opacity']= float(round(array[array<1].max(), n_rounding_places))
    summary_dict['min_opacity']= float(round(array.min(), n_rounding_places))
    summary_dict['stdev_opacity']= float(round(array.std(), n_rounding_places))
    summary_dict['average_opacity']= float(round(array.mean(), n_rounding_places))
    opacity_sum = float(round(array.sum(), n_rounding_places))
    summary_dict['opacity_sum']= float(round(opacity_sum, n_rounding_places))
    
    #Calculate other metrics
    summary_dict['opaque_region_sum']= float(round(array[array<1].sum(), n_rounding_places))
    array_pixel_count = float(round((array.shape[0]*array.shape[1])))
    summary_dict['opaque_fraction']= float(round((array[array<1].shape[0]) / (array_pixel_count), n_rounding_places))
    summary_dict['transmittance']= float(round(opacity_sum / array_pixel_count, n_rounding_places))
    summary_dict['absorbance']= float(round(-1 * math.log10(opacity_sum / array_pixel_count), n_rounding_places+2))
    
    #Generate histogram values and make into JSON-writeable lists
    cnts, bins = np.histogram(array, bins = n_histogram_bins)
    cnts = cnts.tolist()
    bins = np.round(bins[1::], n_rounding_places).tolist()
    summary_dict['histogram_bins']= bins
    summary_dict['histogram_cnts']= cnts
    
    return summary_dict