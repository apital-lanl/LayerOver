# -*- coding: utf-8 -*-
"""
Created:   2024-09-12
Modified:  2025-04-17
Version:   1.1.0 

@author: Aaron Pital (Los Alamos National Lab)

Description: Combined library for gcode, STL, voxelization, and related operations for Direct Ink Write (DIW) work.

"""
version = '1.1.0'
last_modified_date = '2025-04-17'

  #System and built-ins
import os
# import sys
# import json
import math
# import csv
from tkinter import Tk, filedialog

  #Visualizaiton
import matplotlib.pyplot as plt
# from matplotlib import cm, colors
from matplotlib import animation
# from mpl_toolkits.mplot3d import Axes3D

  #Data Handling
import numpy as np
# import pandas as pd

  #Scientifiic algorithm packages
# import scipy.interpolate
# from scipy.spatial import KDTree
# from scipy.spatial import ConvexHull
# from scipy.spatial import Delaunay

  #Get random biz
# from random import choice
# import random

  #Utilities
from tqdm.auto import tqdm

  #LayerOver imports
from Gcode import Gcode
from Camera import Camera
from Point_Clod import Point_Clod
from Opacity import gneral 

# Define package-variables
color_lists = {
'Rainbow':['salmon','darkorange','gold','limegreen',\
           'mediumturquoise','slateblue','mediumorchid'],
'rRainbow': ['salmon','darkorange','gold','limegreen',\
           'mediumturquoise','slateblue','mediumorchid'].reverse(),
}
    

#################################################################################################################
def opacity_from_gcode(filenames_lists = [], n_voxel_points = 5, n_pixels = 100, camera_radius = 1,
                    threshold = 1, max_opacity= 0.1, opacity_func = 'fixed-radial'):
    
    ''' v1.0.0  created: 2024-10-06  modified:  2025-07-02
    
    Description: Tale layer gcode files and a configuration file (LayerUp style)
    
    INPUT:
          #voxel/camera parameters
        n_voxel_points = Number of points for that get evaluated for opacity in a part
        n_pixels =       Number of pixels in 'Camera'; Number of pixels in final opacity measure
        camera_radius =  Radius of points around voxel center to consider
          #Opacity measure parameters
        threshold =      Distance between pixel vector and layer point (in mm typically) 
        max_opacity =    Max decrement multiplier of light crossing through thickest part of strand
        opacity_func =   Function to use for calculating opacity from distance-to-strand-center normal measure
            'fixed-radial'- assumes a constant, round shape 
    '''
    
    
    
    if len(filenames_lists) < 1:
        continue_check = True
        filename_lists = []
        while continue_check:
            # Get the names by user selection
            root = Tk()
            these_filenames = filedialog.askopenfilenames(title= 'Select images to stitch')
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
    if opacity_func == 'fixed-radial':
        opacity_func = lambda distance: max_opacity * (math.sqrt(strand_radius**2 - distance**2) / strand_radius)
    #Default to 'fixed-radial'
    else:
        opacity_func = lambda distance: max_opacity * (math.sqrt(strand_radius**2 - distance**2) / strand_radius)
    
    #%% (3) Do everything
    for file_list in filename_lists:
        try:
            part = Gcode.Gcode()
            part.get_gcode(files = file_list)
            
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
                
            #################################################################################################################
            #Generate center point, normals, etc. from fibonacci points and near neighbors
            point_dicts = Camera.get_random_point_summaries(part, n_retrieval_points = n_voxel_points, \
                                               excluded_point_indcs = [], show_iterations = True)
                
            #################################################################################################################
            # Make voxels and visualize
            
            voxel_opacity_images = []
            for idx, key in enumerate(list(point_dicts.keys())):
                
                #generate voxel name
                voxel_name = str(part_name_guess + "_voxel-" + str(idx))
                
                #################################################################################################################
                #pull voxel features
                center = point_dicts[idx]['approximate_center']
                normal_vector = point_dicts[idx]['unit_normal']
                radius = camera_radius      #TODO: leave option to make this layer-defined
                #TODO: Radius is currently hard-coded by user input (Cell 1); could be made to be more flexible here
                
                #get voxel bounds and add to dictionary 
                voxel_dict = Camera.get_voxel(center, normal_vector, radius= radius, voxel_depth= 2*radius, n_pixels=n_pixels)
                voxel_dict = Camera.get_layer_points(part, voxel_dict, boundary_buffer_multiplier = 0.2)
                
                #################################################################################################################
                #PLot voxel
                gif = False
                gif_frames = 4     #number of frames in GIF
                gif_integrals = 50   #delay between frames in msec
                azim_total_degrees = 360
                azim_step = azim_total_degrees/gif_frames
                
                plot_first_layer = True
                plot_last_layer = True
                
                top_pixel_coordinate_array = voxel_dict['top_pixel_coordinate_array']
                bottom_pixel_coordinate_array = voxel_dict['bottom_pixel_coordinate_array']
                
                bottom_layer_name = layer_list[0]
                top_layer_name = layer_list[-1]
                
                first_points = np.array([0])
                for c_idx, contig_list in enumerate(voxel_dict['layer_coordinates'][bottom_layer_name]['coordinates']):
                    for p_idx, point in enumerate(contig_list):
                        this_point = point
                        #ignore 0-dim arrays
                        if this_point.shape[0]<3:
                            pass
                        elif len(first_points.shape)<2:
                            first_points = np.array([this_point])
                        else:
                            first_points = np.append(first_points, [this_point], axis =0)
                first_points = np.transpose(first_points)
                
                last_points = np.array([0])
                for c_idx, contig_list in enumerate(voxel_dict['layer_coordinates'][top_layer_name]['coordinates']):
                    for p_idx, point in enumerate(contig_list):
                        this_point = point
                        #ignore 0-dim arrays
                        if this_point.shape[0]<3:
                            pass
                        elif len(last_points.shape)<2:
                            last_points = np.array([this_point])
                        else:
                            last_points = np.append(last_points, [this_point], axis =0)
                last_points = np.transpose(last_points)
                
                  #split x/y/z for plotting
                top_xs = top_pixel_coordinate_array[::, 0]
                top_ys = top_pixel_coordinate_array[::, 1]
                top_zs = top_pixel_coordinate_array[::, 2] 
                bot_xs = bottom_pixel_coordinate_array[::, 0]
                bot_ys = bottom_pixel_coordinate_array[::, 1]
                bot_zs = bottom_pixel_coordinate_array[::, 2] 
                
                #Initialize 3D plot
                fig = plt.figure(figsize=(10,10))
                ax = fig.add_subplot(111, projection = '3d')
                
                  #add substrate points
                ax.scatter(top_xs, top_ys, top_zs, s=.2, color = 'teal', alpha=0.5)
                ax.scatter(bot_xs, bot_ys, bot_zs, s=.2, color = 'brown', alpha=0.5)
                
                if plot_first_layer:
                    ax.scatter(first_points[0], first_points[1], first_points[2], s=10, color = 'blue', alpha=0.1)
                if plot_last_layer:
                    ax.scatter(last_points[0], last_points[1], last_points[2], s=10, color = 'purple', alpha=0.1)
                
                def init():
                    ax.view_init(elev=10., azim=0)
                    return [fig]
                
                def animate(i):
                    global azim_step
                    ax.view_init(elev=10., azim=i*azim_step)
                    print(f"Running {i+1} of {gif_frames}")
                    return [fig]
                
                # Animate
                anim = animation.FuncAnimation(fig, animate, init_func=init,
                                               frames=gif_frames, interval=gif_integrals, blit=True)
                if gif:
                    # Save; Takes a while (~5 minutes?)
                    save_format = ".mp4"
                    format_name = save_format.replace('.','').upper()
                    
                    root = Tk()
                    save_filename = filedialog.asksaveasfilename(filetypes=[(format_name, save_format)])
                    root.destroy()
                    
                    anim.save(save_filename+save_format, fps=30, extra_args=['-vcodec', 'libx264'])
                    #anim.save(save_filename, fps=30)
                else:
                    plt.title(voxel_name)
                    voxel_point_name = str(part_name_guess + "_VoxelPoints" + str(idx))
                    save_name = os.path.join(part_dir, voxel_point_name)
                    plt.savefig(save_name, dpi=600)
                    plt.show()
                    
                #################################################################################################################
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
                
                #################################################################################################################
                #Get the actual distances for each pixel and calculate opacity
                pixel_min_distances = []
                opacity_array = np.ones( (voxel_dict['top_pixel_coordinate_array'].shape[0]))
                for l_idx, layer_key in enumerate(list(part.layers.keys())):
                
                    flat_coordinates = voxel_dict['flat_layer_coordinates'][layer_key]['coordinates']
                    flat_indices = voxel_dict['flat_layer_coordinates'][layer_key]['indices']
                
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
                        distances = np.array([Point_Clod.line_line_distance(layerp1, layerp2, bot_pixel, top_pixel)[2]
                                     for layerp1, layerp2 in close_segment_points])
                        min_distance = np.min(distances)
                        
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
                        opacity_array[pixel_idx] = opacity_array[pixel_idx] * this_opacity
                    
                opacity_image = np.reshape(opacity_array, (n_pixels, n_pixels))
                voxel_opacity_images.append(opacity_image)
                
                plt.imshow(opacity_image)
                plt.title(voxel_name)
                save_name = os.path.join(part_dir, str(voxel_name+ '.png'))
                plt.imsave(save_name, opacity_image, dpi=600)
                plt.show()
            
            opacity_arrays = np.array(voxel_opacity_images)
            partial_name = str(part_name_guess + '_OpacityArrays.npy')
            save_name = os.path.join(part_dir, partial_name)
            np.save(save_name, opacity_arrays)
        except IndexError as IE:
            print(IE)
            print()
            print("Random failure")    
        
        
###############################################################################################################################       
    

    
    
###############################################################################################################################
    
    

###############################################################################################################################

#%%

