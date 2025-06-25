# -*- coding: utf-8 -*-
"""
Created:   2024-04-10
Modified:  2024-05-12
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Object class for POV renders, POV plots, and geometry-based measures of a collection of XYZ cartesian geometries.

"""
version = '0.0.0'
last_modified_date = '2000-00-00'

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

  #Other LayerOver modules
from Point_Clod import Point_Clod


class Camera():

    def __init__(self, gcode_part, lens_dict = None, render_dict= None, view_direction = 'in'):
        ''' v0.2.0   created:2024-04-10   modified:2024-05-08

        INPUTS:  'cart_location_dict' - dictionary of points to view from; need location and normal vector
                 'render_dict' - contains all the geometry information
                 'strand_geometry' - 'by-layer'-- each layer has it's own strand geometry
        '''

        self.gcode_part = gcode_part
        
        # If 'lens_dict' isn't provided, make one
        if lens_dict == None:
            self.lens_dict = {
                'type': 'simple',
                'shape': 'square',
                'profile': 'flat',
                'size': (6,6),
                'working_distance': 6
                }
        else:
            self.lens_dict = lens_dict

          #Modify with gcode_part
        

        # If 'render_dict' isn't specified, make one
        if render_dict == None:
            self.render_dict = {
                'geometry_dict':
                    {
                    'strand_size_type': 'global',
                    'strand_diameter': 150,
                    },
                
                'intensity_cmap': None,
                'perspective_padding': 1
                }

        #Instantiate view dictionaries
        self.povs = {}
        self.voxels = {}

    
    def make_povs(self, number_of_locations, location_list =None):
        ''' v0.1.0   created:2024-05-09   modified:2024-05-20

        INPUTS:  'location_number' - number of points to pull cameras on
                 (optional)
                 'location_list'- supply a list of discrete locations rather than pick at random
        '''

        if location_list == None:
        
            point_dicts = Camera.get_random_point_summaries(self.gcode_part, n_retrieval_points = number_of_locations, \
                                                            excluded_point_indcs = [], show_iterations = True)
            self.pov_point_summaries = point_dicts
        else:
            #TODO: make this part work
            pass

        for idx, key in enumerate(list(point_dicts.keys())):
            center = point_dicts[idx]['approximate_center']
            normal_vector = point_dicts[idx]['unit_normal']
            radius = self.lens_dict['size'][0]
            
            #get voxel bounds and add to dictionary 
            voxel_dict = Camera.get_voxel(center, normal_vector, radius= radius, voxel_depth= 2*radius, n_pixels=1000)
            voxel_dict = Camera.get_layer_points(self, voxel_dict)

            #Updtate 'voxel' dictionary on object
            self.voxels.update( {key:voxel_dict} )
    
    
    def populate_povs(self, make_opacity = True, direction = 'in'):
        
        ''' v0.1.1   created:2024-05-09   modified:2024-05-11

        INPUTS:  'cart_location_dict' - dictionary of points to view from; need location and normal vector
                 'render_dict' - contains all the geometry information
                 'strand_geometry' - 'by-layer'-- each layer has it's own strand geometry
        '''
        
        #Initialize the biz
        point_id_list = list(self.location_dict.keys())
        
        #Run through each point's properties in 'cart_location_dict', get properties and apply to 'self.povs' dictionary
        for point_id in point_id_list:
            this_center = self.location_dict[point_id]['approximate_center']
            this_normal = self.location_dict[point_id]['unit_normal']
            #Derive eqn of plane (cx*nx + cy*ny + cz*nz = d) by plugging in center to normal
            plane_d = this_center[0]*this_normal[0] + this_center[1]*this_normal[1] + this_center[2]*this_normal[2]
              
              #Create
            if direction == 'in':
                pass
            elif direction == 'out':
                pass
            else:
                pass    #default to 
            
            
            if self.render_dict['geometry_dict']['strand_size_type'] == 'global':
                pass


    


    ############################################################################################################################
    ####  Utilities and StaticMethods  #########################################################################################
    ############################################################################################################################

    @staticmethod
    def get_random_point_summaries(gcode_part, n_retrieval_points = 5, \
                                   excluded_point_indcs = [], show_iterations = False):
    
        ''' v0.2.2  created:2024-05-08  modified:2024-05-12 
        
        INPUT:   'gcode_part'- a part generated from 'gcode' class; inclueds config.JSON and .PGM layer files
            (optional)
                 'n_retrieval_points' - number of points to evaluate and return
                 
        ACTION:  Pull fibonacci point triangles from part and calculate center, normal, and other stats 
        OUTPUT:  Dictionary with each point's summary data; dictionary keys are just integer #'s starting at 0 (i.e. index of point_
        '''
        
        #v0.2.0 version; takes known fibonacci points (that seem like they work) and makes a quick 'nearest neighbor' triangle to other fib. points
          #Get stats about part and initialize variables
        point_n_total = gcode_part.substrate_dict['fib_points'].shape[0]
        eval_point_dict = {}
        
        #Get 'n' triangles from the surface randomly
        rand_point_indcs = []
        rand_triang_coord = []
        for i in range(0, n_retrieval_points):
            this_choice = choice(range(0, point_n_total))
            this_point = gcode_part.substrate_dict['fib_points'][this_choice]
              #add to list if it's not already there (don't want duplicates)
            these_coords = [this_point]
            if (this_choice not in rand_point_indcs) and (this_choice not in excluded_point_indcs):
                rand_point_indcs.append(this_choice)
                trial_points=[]
                for idx, trial_point in enumerate(gcode_part.substrate_dict['fib_points']):
                    this_distance = math.sqrt((trial_point[0]-this_point[0])**2 + (trial_point[1]-this_point[1])**2 + (trial_point[2]-this_point[2])**2)
                    if this_distance != 0:
                        trial_points.append((idx, this_distance))
                #Sort coordinates by distance and add two nearest points
                sorted_trial_points = sorted(trial_points, key=lambda x: x[1])
                  #smallest values sorted first; pick the smallest two
                second_point = gcode_part.substrate_dict['fib_points'][sorted_trial_points[0][0]]
                these_coords.append(second_point)
                  #calculate angle betwen p0->p1 and p0->p2; if about 180 (or 0), use a new p2
                third_point = gcode_part.substrate_dict['fib_points'][sorted_trial_points[1][0]]
                vector_1 = second_point-this_point
                vector_2 = third_point-this_point
                this_angle = math.acos( ((vector_1[0]*vector_2[0]+vector_1[1]*vector_2[1]+vector_1[2]*vector_2[2])/  \
                            (math.sqrt(vector_1[0]**2+vector_1[1]**2+vector_1[2]**2) * \
                             math.sqrt(vector_2[0]**2+vector_2[1]**2+vector_2[2]**2))))
                if this_angle >2:
                    #try another point; shouldn't be possible to be co-linearish at this point
                    third_point = gcode_part.substrate_dict['fib_points'][sorted_trial_points[2][0]]
                    vector_1 = second_point-this_point
                    vector_2 = third_point-this_point
                    this_angle = math.acos( ((vector_1[0]*vector_2[0]+vector_1[1]*vector_2[1]+vector_1[2]*vector_2[2])/  \
                            (math.sqrt(vector_1[0]**2+vector_1[1]**2+vector_1[2]**2) * \
                             math.sqrt(vector_2[0]**2+vector_2[1]**2+vector_2[2]**2))))
                    these_coords.append(third_point)
                else:
                    these_coords.append(third_point)
                  #Add triangle coordinates to list
                rand_triang_coord.append(these_coords)
        
        ''' v0.1 version; uses convex hull function that I don't trust fully
        #Get stats about part and initialize variables
        triangle_n_total = gcode_part.substrate_dict['fib_surface_triangles'].shape[0]
        eval_point_dict = {}
        
        #Get 'n' triangles from the surface randomly
        rand_triang_indcs = []
        rand_triang_points = []
        rand_triang_coord = []
        for i in range(0, n_retrieval_points):
            this_choice = choice(range(0, triangle_n_total))
            these_indcs = gcode_part.substrate_dict['fib_surface_triangles'][this_choice]
              #add to list if it's not already there (don't want duplicates)
            these_coords = []
            if (this_choice not in rand_triang_indcs) and (this_choice not in excluded_point_indcs):
                rand_triang_indcs.append(this_choice)
                rand_triang_points.append(these_indcs)
                for idx in these_indcs:
                    these_coords.append(gcode_part.substrate_dict['fib_points'][idx])
                rand_triang_coord.append(these_coords)
        '''
        
        #Run through each set of triangles and get normal, center, equation of plane, etc. for them
        point_dicts = {}
        for idx, point_coords in enumerate(rand_triang_coord):
            this_point_dict = Point_Clod.get_3point_normal(point_coords, show_plot = show_iterations)
            point_dicts.update({idx:this_point_dict})
    
        return point_dicts


    @staticmethod
    def get_layer_points(gcode_part, voxel_dict, boundary_buffer_multiplier = 0.1):
        
        ''' v0.2.0  created:2024-05-09  modified:2024-06-03
    
        Take a part with layers and a dictionary of voxel properties and return a dictionary of layer-points within
          the voxel region and some properties.
        
            INPUT:   'gcode_part'- gcode class instance with layers added (happens at initialization if 'Layer -' files are selected)
                     'voxel_dict'- dictionary with voxel bounds, etc.
            ACTION:  cull points to within voxel bounds; easy because they're square at this point (that's a joke)
            OUTPUT:  updated 'voxel_dict' with layers and in-bound coordinates added as lists;
                        NOTE: lists have inhomogenous dimensions, so naive parsing is kind of a pain.
        
        '''
    
        #Initialize new entry for layer coordinates
        voxel_dict.update( {'layer_coordinates':{}})
        
        #TODO: make sure the output of the 'voxel_bounding_box' is ordered somehow
        bounds_array = []
        for key in list(voxel_dict['voxel_bounding_box']):
            bounds_array.append(voxel_dict['voxel_bounding_box'][key])
        xmin, xmax, ymin, ymax, zmin, zmax = bounds_array
        xmin -= abs(xmax-xmin)*boundary_buffer_multiplier
        xmax += abs(xmax-xmin)*boundary_buffer_multiplier
        ymin -= abs(ymax-ymin)*boundary_buffer_multiplier
        ymax += abs(ymax-ymin)*boundary_buffer_multiplier
        zmin -= abs(zmax-zmin)*boundary_buffer_multiplier
        zmax += abs(zmax-zmin)*boundary_buffer_multiplier
    
        #run through each layer's coordinates and cull to dictionary list
        for layer_name in list(gcode_part.layers.keys()):
              #initialize dictionary entry for this layer
            voxel_dict['layer_coordinates'].update( {layer_name:{}} )
            
            this_layer_dict = gcode_part.layers[layer_name]
            these_coordinates = this_layer_dict['coordinates']
            #save as lists for appending to dictionary
            in_bound_indcs = []
            in_bound_coordinates = []
            seq_coordinates = []
            last_idx_hit = -1  #placeholder for denoting sequential coordinates
            for idx, coordinates in enumerate(these_coordinates):
                x, y, z = coordinates
                
                if (x>xmin) and (x<xmax) and (y>ymin) and (y<ymax) and (z>zmin) and (z<zmax):
                    in_bound_indcs.append(idx)
                    #want to keep sequential coordinates nested
                    if (idx-last_idx_hit) > 1:
                        in_bound_coordinates.append(seq_coordinates)  #add the previous list
                        seq_coordinates = [coordinates]  #start a new list with this index's coordinates
                    else:
                        seq_coordinates.append(coordinates)
                        
                    last_idx_hit = idx
            in_bound_coordinates.append(seq_coordinates)  #append whatever is left in the queue
    
            #append lists to dictionary under 'layer_name' for later pulling from gcode
            voxel_dict['layer_coordinates'][layer_name].update( {'indices': in_bound_indcs} )
            voxel_dict['layer_coordinates'][layer_name].update( {'coordinates': in_bound_coordinates} )
                
        return voxel_dict
    
    
    @staticmethod
    def make_camera(center_point, unit_normal, \
                    plane_size = (6,6), lens_type = 'square', \
                    focal_pyramid_length = 6, \
                    camera_color= 'k', opacity = 0.7 \
                    ):
        ''' v0.1.0  created:2024-05-10  modified:2024-05-10

        Take 
        
            INPUT:   'center'- point at center of camera/lens plane
                     'unit_normal'- unit normal vector for direction of camera
            ACTION:  math
            OUTPUT:  'camera_dict'-
                        'plane_points'- 4 pts (x,y,z) describing render plane
                        'lens_type'- defaults to a square lens 6mm x 6mm
                        'voxel_points'- 8 pts (x,y,z) describing bounding box for camera's observation voxel
                             (includes padding so geometries behave)
                        'focal_pyramid_point'- point 'behind' camera lens that defines focal cone (shown as pyramid)
        '''
        
        #plane
        
        #return camera_dict
        
        pass

    
    @staticmethod
    def get_voxel(center, vector, radius= 3, voxel_depth= 6, n_pixels=100):
        
        ''' v0.2.4  created:2024-05-09  modified:20245-04-17
        v0.2.3- modified bounding points to be actual 'pixel' array rather than radial points
        v0.2.4- added 'top' or 'bottom' voxel face check
        
        Given a point and a vector, generate two circles at two points along with 'center' in middle.
            Generate a set of min:max bounds defining rectangular voxel for that cylindrical region.
            Will be bounding box for region, not true voxel region.
            INPUT:   'vector'- should be normal vector to plane
            ACTION:  lorem
            OUTPUT:  lorem
            
        Changes:
            2025-04-17- changed 'n_circ_points' to reflect change in function. Used to be (360/n_circ_points)= "num of points returned"; now, n_circ_points is simply the "num of points returned"
        '''
    
        #Norm vector just in case it's not already
        vector /= np.linalg.norm(vector)
        
        #Get two bounding points for cylindrical region (center of each end's face); 'top'/'bottom' are arbitrary
        top_center = center + (vector*(voxel_depth/2))
        bottom_center = center - (vector*(voxel_depth/2))
          #Generate points for circumerence of cirlce
          #   will be 360/'n_circ_points' number of points
        top_radial_points = Point_Clod.generate_radial_points(vector, top_center, radius, n_circ_points = 6)
        bottom_radial_points = Point_Clod.generate_radial_points(vector, bottom_center, radius, n_circ_points = 6)
    
        #Get corners of visualization square
        pixel_array_radius = math.sqrt(2*radius**2)
        top_pixel_array_corners = Point_Clod.generate_radial_points(vector, top_center, pixel_array_radius, n_circ_points = 4)
        bottom_pixel_array_corners = Point_Clod.generate_radial_points(vector, bottom_center, pixel_array_radius, n_circ_points = 4)
        
        #Draw 'top' pixel array by coordinates
        pa = top_pixel_array_corners[0]
        pb = top_pixel_array_corners[1]
        pc = top_pixel_array_corners[2]

        # Generate point-wise grid based on plane equation; 
        pixel_coordinate_array = np.zeros((3, n_pixels**2))
        idx = 0
        for u in tqdm(range(0, n_pixels)):
            for v in range(0,n_pixels):
                this_point = np.multiply((u/n_pixels), pa) + \
                             np.multiply((v/n_pixels), pb) + \
                             np.multiply(((n_pixels-u-v)/n_pixels), pc)
                x,y,z = this_point
                pixel_coordinate_array[0,idx] = x
                pixel_coordinate_array[1,idx] = y
                pixel_coordinate_array[2,idx] = z
                idx += 1
        temp_top_pixel_coordinate_array = pixel_coordinate_array.transpose()
        
        #Draw 'bottom' pixel array by coordinates
        pa = bottom_pixel_array_corners[0]
        pb = bottom_pixel_array_corners[1]
        pc = bottom_pixel_array_corners[2]

        # Generate point-wise grid based on plane equation; 
        pixel_coordinate_array = np.zeros((3, n_pixels**2))
        idx = 0
        for u in tqdm(range(0, n_pixels)):
            for v in range(0,n_pixels):
                this_point = np.multiply((u/n_pixels), pa) + \
                             np.multiply((v/n_pixels), pb) + \
                             np.multiply(((n_pixels-u-v)/n_pixels), pc)
                x,y,z = this_point
                pixel_coordinate_array[0,idx] = x
                pixel_coordinate_array[1,idx] = y
                pixel_coordinate_array[2,idx] = z
                idx += 1
        temp_bottom_pixel_coordinate_array = pixel_coordinate_array.transpose()

        #Check which array is 'outside' (further away from origin) and make that 'top_array'
          #pick first points and measure distance from origin
        ptopx, ptopy, ptopz = temp_top_pixel_coordinate_array[0]
        pbotx, pboty, pbotz = temp_bottom_pixel_coordinate_array[0]
        top_distance = math.sqrt(ptopx**2 + ptopy**2 + ptopz**2)
        bot_distance = math.sqrt(pbotx**2 + pboty**2 + pbotz**2)
        if bot_distance > top_distance:
            bottom_pixel_coordinate_array = temp_top_pixel_coordinate_array
            top_pixel_coordinate_array = temp_bottom_pixel_coordinate_array
            print("Flipped: ")
            print(f"Top: {temp_bottom_pixel_coordinate_array[0]}  {bot_distance}")
            print(f"Bot: {temp_top_pixel_coordinate_array[0]}  {top_distance}")
            print()
        else:
            bottom_pixel_coordinate_array = temp_bottom_pixel_coordinate_array
            top_pixel_coordinate_array = temp_top_pixel_coordinate_array
            print("Not Flipped: ")
            print(f"Top: {temp_top_pixel_coordinate_array[0]}  {top_distance}")
            print(f"Bot: {temp_bottom_pixel_coordinate_array[0]}  {bot_distance}")
            print()
            
        #Whip through circumferential points and get min:max bounding box boundaries 
        x_min = 10000
        x_max = -10000
        y_min = 10000
        y_max = -10000
        z_min = 10000
        z_max = -10000
        
        #used to be 'top_radial_points'
        for point in top_pixel_coordinate_array:
            if point[0] < x_min:
                x_min = point[0]
            if point[0] > x_max:
                x_max = point[0]
            if point[1] < y_min:
                y_min = point[1]
            if point[1] > y_max:
                y_max = point[1]
            if point[2] < z_min:
                z_min = point[2]
            if point[2] > z_max:
                z_max = point[2]
        for point in bottom_pixel_coordinate_array:
            if point[0] < x_min:
                x_min = point[0]
            if point[0] > x_max:
                x_max = point[0]
            if point[1] < y_min:
                y_min = point[1]
            if point[1] > y_max:
                y_max = point[1]
            if point[2] < z_min:
                z_min = point[2]
            if point[2] > z_max:
                z_max = point[2]
    
        #Make a return_dict
        voxel_dict = {
            'top_cylinder_circum_points': top_radial_points,
            'bottom_cylinder_circum_points': bottom_radial_points,
            'top_pixel_array_corners':  top_pixel_array_corners,
            'bottom_pixel_array_corners':  bottom_pixel_array_corners,
            'top_pixel_coordinate_array': top_pixel_coordinate_array,
            'bottom_pixel_coordinate_array': bottom_pixel_coordinate_array,
            #TODO: add plane equations
            #point_clod.get_3point_normal(point_coords)
            'voxel_bounding_box': {
                'xmin': x_min,
                'xmax': x_max,
                'ymin': y_min,
                'ymax': y_max,
                'zmin': z_min,
                'zmax': z_max},
            'center': center,
            'vector': vector,
            'voxel_vertices': {
                 1: [x_min, y_min, z_min],
                 2: [x_max, y_max, z_min],
                 3: [x_min, y_max, z_min],
                 4: [x_max, y_min, z_min],
                 5: [x_min, y_min, z_max],
                 6: [x_max, y_max, z_max],
                 7: [x_min, y_max, z_max],
                 8: [x_max, y_min, z_max]
                    },
            'voxel_faces': {
                'bottom':[1,3,2,4],
                'left':[1,3,7,5],
                'top':[5,8,6,7],
                'front':[1,4,8,5],
                'right':[4,2,6,8],
                'back':[2,6,7,3]
                }
                    }
    
        return voxel_dict


