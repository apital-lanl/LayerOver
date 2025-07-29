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
import traceback

  #Visualizaiton
import matplotlib.pyplot as plt

  #Data Handling
import numpy as np

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
    'max_opacity': 0.1,             # max decrement multiplier of light crossing through thickest part of strand
    }

general_opacity_func = lambda distance: max_opacity * (math.sqrt(strand_radius**2 - distance**2) / strand_radius)



def opacity_from_gcode(filenames_lists = [], n_voxel_points = 5, n_pixels = 100, camera_radius = 3,
                    threshold = 1, max_opacity= 0.1, opacity_func = 'fixed-radial', size_multiplier = 1,
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

    squish_factor= 0.3          #compaction ratio vs. ideal stack (0.3 = 30% compaction, i.e. 1" part ends up being 0.7")
    base_compaction_ratio= 0.2  #ratio of total 'squish_factor' that occurs at the baseplate; i.e. '0.3' means 30% of normal compaction occurs at baseplate
    
    if len(filenames_lists) < 1:
        continue_check = True
        filename_lists = []
        while continue_check:
            # Get the names by user selection
            root = Tk()
            these_filenames = filedialog.askopenfilenames(title= 'Select images to stitch', filetypes = [('LayerOver Files', '*.json *.pgm')])
            root.destroy()
            
            # Sort to get images in left-to-right, top-to-bottom order by index (hopefully)
            filename_lists.append(these_filenames)
            
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
    for file_list in filename_lists:
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
                        if homogenous_strand_size:
                            pass  #strand_radius will already be defined
                        else:
                            strand_radius = radius_list[l_idx]
                        
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

                    partial_layer_name = str(part_name_guess + f"_{voxel_name}-{layer_key}.npy")
                    save_name = os.path.join(part_dir, partial_layer_name)
                    np.save(save_name, layer_opacity_image)
                
                #Combine each layer opacity into a single voxel's opacity
                #  Adjusts for the fact that strands coallesce to some extent
                for layer_idx, layer_image in enumerate(layer_opacity_images):
                    if layer_idx == 1:
                        max_squished_opacity = max_opacity * (1-(squish_factor*base_compaction_ratio))  #assume baselayer squishes all over at max-strand-thickness, by about 'base_compaction_ratio'
                        layer_image[layer_image>max_squished_opacity]= max_squished_opacity
                        opacity_array = layer_image
                    elif layer_idx ==2:
                        sum_image = layer_opacity_images[layer_idx-1] * layer_image
                        max_squished_opacity = max_opacity**2 * (1-squish_factor)  #Each overlapping region will be (1-squish_factor) less opaque than expected 

                    else:


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

    voxel_dict.update({'point_dicts':point_dicts})
    voxel_dict.update({'radius_list': radius_list})

    return voxel_dict


def opacity_from_voxel_dict(voxel_dict = None, file_list = None, n_voxel_points = 1, n_pixels = 100, camera_radius = 3,
                    threshold = 1, max_opacity= 0.1, opacity_func = 'fixed-radial', size_multiplier = 1,
                    strand_thickness_func = 'circular'):
    
    ''' v0.1.0  created: 2025-07-28  modified:  2025-07-28
    
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


def adjust_opacity_for_squish(layer_array, combined_opacity_array,
                              full_strand_opacity = 0.1,
                              opacity_function = 'fixed-radial'):





    return new_combined_opacity_array