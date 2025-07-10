# -*- coding: utf-8 -*-
"""
Created:   2024-03-18
Modified:  2024-03-22
Version:   0.4.0

@author: Aaron Pital (Los Alamos National Lab)

Description: 

"""
version = '0.4.0'
last_modified_date = '2024-03-22'

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

  #Import other LayerUp modules
from Point_Clod import Point_Clod

  #Utilities
from tqdm.auto import tqdm

class Gcode():

    '''  
    Summary: 
    ''' 
    
    def __init__(self, config=None):
        ''' v1.1  created:2024-03-17  modified:2024-03-22
        '''
        # Initialize Class with some generic placeholders to be filled in later
        self.layer_props = {
            '1': {}
        }
        self.spatial_params = ['X', 'Y', 'Z']
        self.function_params = ['DISPA']
        self.offsets ={}
        self.registration_points = {}
        self.layers = {}
        self.dir_path = ''
        self.config_path = ''
        self.layer_files = []
        self.substrate_dict = {}

    
    def get_gcode(self, files= '', add_layers=True):
        ''' v1.0  created:2024-03-17  modified:2024-03-20
        Calculate Euclidean distance and cull points that are outside tolerance range
        '''
        if files == '':
            # Get LayerUp files; JSON files are config files and PGM are the Gcode for each layer
            root = Tk()
            files = filedialog.askopenfilenames(title="Select gcode and config files", \
                                                filetypes=[("LayerUp Files",'*.json *.pgm')])
            root.destroy()

        # Go through config files and check for uniqueness
        config_guesses = [file for file in files if 'config' in os.path.basename(file).lower()]
        if len(config_guesses) == 1:
            self.config_path = config_guesses[0]
            self.parse_config(self.config_path)
        else:
            print(config_guesses)
        
        # Collect Gcode
        for file in files:
            if file.endswith(".pgm") or file.endswith(".txt"):
                self.layer_files.append(file)

        if add_layers:
            for layer_filename in self.layer_files:
                layer_name = os.path.basename(layer_filename.replace('.pgm',''))
                xyzacf_array = Gcode.parse_gcode_file(layer_filename, part = self)
                layer_dict = Point_Clod.AC_to_XYZ(xyzacf_array, self.substrate_dict, offset_adjust=False)
                layer_dict.update( {'layer_name': layer_name} )
                self.layers.update( {layer_name: layer_dict} )
                print(f"Layer added: {layer_name}")
                print()

        
    def parse_config(self, filepath):
        ''' v1.0  created:2024-03-17  modified:2024-03-20
        Get values and parameters from 
        '''
        with open(filepath, 'r', encoding ='utf-8') as file:
            config_dict = json.load(file)

        self.config_dict = config_dict
        
        #Get strand diameter from config data
        layer_number = len(config_dict['part']['layers']['value'])
        diameter_list = []
        unique_diameters = []
          #pull data from the structure of each layer
          #TODO: only the strand diameter is pulled from this, but pitch, polar angle, etc. could be pulled as well
        for idx in range(layer_number):
            this_diameter = config_dict['part']['layers']['value'][idx]['structure']['value']['strand_diameter']['value']
            diameter_list.append(this_diameter)
            if this_diameter not in unique_diameters:
                unique_diameters.append(this_diameter)
        
        self.layer_strand_diameters = diameter_list
        if len(unique_diameters)==1:
            self.strand_diameter = unique_diameters[0]
        else:
            self.strand_diameter = 0

        #'Config' substrate parsing as of LayerUp 2.1.1 for KC 5-axis machine
        substrate_flag_guess = config_dict['substrate']['classname']
          #Parse 'config_dict' data to validate guess at substrate
        if substrate_flag_guess == 'Mandrel' or substrate_flag_guess =='FlatPlate':
            substrate_flag = substrate_flag_guess
        
          #For flat plate
        if substrate_flag == 'FlatPlate':
              #Make generic X-Y plate by dimensions
            length = config_dict['substrate']['length']['value']
            width = config_dict['substrate']['width']['value']
              #Add fibonacci points from center of plate; clip at edges (strange edge-point density)
            sub_dict = self.make_flatplate_xyzpoints(width, length)
        
          #For 'parameter'-driven mandrel shape
        if substrate_flag == 'Mandrel':
              #Pull shape profile from 'lats' and 'verts' (equivalent to X and Y)
            lats = config_dict['substrate']['shapeFunction']['value']['lat']
            verts = config_dict['substrate']['shapeFunction']['value']['vert']
              #Take mandrel shape and project fibonacci points onto it
            sub_dict = self.make_mandrel_xyzpoints(lats, verts)
              # Assume a mandrel is on a 5-axis and add the Azimuthal and rotational(C) parameters
            self.spatial_params.append('A')
            self.spatial_params.append('C')
        
        self.substrate_dict = sub_dict

    
    ############################################################################################################################
    ####  Utilities and StaticMethods  #########################################################################################
    ############################################################################################################################
    
    @staticmethod
    def make_flatplate_xyzpoints(width, length, offsets = None, show = True, gif_make=False):
        ''' v1.0  created:2024-03-17  modified:2024-03-21
        Generate a flat plate from width and length from config.json file 
        '''
        # Calculate the minimum circle radius needed to bound the plate
        #   this is just the hypoteneuse of the plate triangle/2 
        plate_hypot = math.sqrt(width**2 + length**2)
        min_radius = plate_hypot/2
        
        lats = np.linspace(0,plate_hypot, 1000)
        if offsets == None:
            verts = lats*0
        else:
            verts = offsets
        
        #lats are X, verts are Z, and Y is interpolated
        xs, zs, ys, radii = Point_Clod.fibonnaci_points(verts, lats)
          #remove points out of the plate bounds; zero is the center
        plate_xs =[]
        plate_ys = []
        plate_zs = []
        for x, y, z in zip(xs, ys, zs):
            if (x< width/2 and x> (-width/2)) and \
               (y< length/2 and y> (-length/2)):
                plate_xs.append(x)
                plate_ys.append(y)
                plate_zs.append(z)

        # Make some triangles
        xyzarray, sphereInds, triangInds = Point_Clod.convexHull(plate_xs,plate_ys,plate_zs)

        # Show the biz
        Point_Clod.visualize_3D(xyzarray, sphereInds, triangInds, gif=gif_make, point_color='k')

        substrate_dict ={
            'width': width,
            'length': length,
            'surface_area': width*length,
            'centerpoint': (0,0),
            'fib_points': xyzarray,
            'fib_surface_triangles': triangInds,
            'fib_sphere_indcs': sphereInds
            }
        
        return substrate_dict

    
    @staticmethod
    def make_mandrel_xyzpoints(lats, verts):
        ''' v1.0  created:2024-03-17  modified:2024-03-21
        Generate a mandrel from config.json file profile points
        '''
        
        xs, ys, zs, radii = Point_Clod.fibonnaci_points(verts, lats)
        
        # Make some triangles
        xyzarray, sphereInds, triangInds = Point_Clod.convexHull(xs,zs,ys)
        
        # plot the test
        Point_Clod.visualize_3D(xyzarray, sphereInds, triangInds, gif=False, point_color ='k') 
        
        substrate_dict = {
            'radius': max(lats)-min(lats),
            'height': max(verts)-min(verts),
            'surface_area': 3*math.pi*(max(lats)-min(lats))**2,
            'centerpoint': (0,0),
            'surface_profile':np.stack((np.linspace(0,90, len(lats)),lats,verts), axis=1),
            'fib_points': xyzarray,
            'fib_surface_triangles': triangInds,
            'fib_sphere_indcs': sphereInds
            }
        return substrate_dict
    
    
    @staticmethod 
    def recursive_dictkey_printer(dict, level):
        '''
        v1.0   created: ?   modified
        '''
        
        keys = list(dict.keys())
        for key in keys:
            tab_space = "   "*level
            if level == 0:
                print()
                print('----------------------------------------------------')
                print()
                print(f"{tab_space}{key}")
            else:
                print(f"{tab_space}{key}")
            try:
                Gcode.recursive_dictkey_printer(dict[key],level+1)
            except:
                #print(dict[key])
                pass
    
    
    @staticmethod
    def parse_gcode_file(filepath, part = ''):

        '''
        v1.0   created: ?   modified
        '''
        
        with open(filepath, 'r', encoding ='utf-8') as file:
            text = file.read().split('\n')
        
        header_flag = False
        toolpath_flag = False
        footer_flag = False
        xyzacf_mode_array = []   #save points to a single list; this will be rolled into a numpy array later, but keep as boring list for now
        idx_dict = {
            'X':  0,
            'Y':  1,
            'Z':  2,
            'A':  3,
            'C':  4,
            'F':  5,
            'G1': 6,
            }
        
        if part=='':
            part = Gcode()
        
        for line in text:
            # Almost nothing starts with ' F', except for feed rate on 5-axis parts
               #add whitespace if X,Y,Z,C is negative; otherwise '-' replaces whtiespace for some reason
            clean_line = line.replace('X-', 'X -')
            clean_line = clean_line.replace('Y-', 'Y -')
            clean_line = clean_line.replace('Z-', 'Z -')
            clean_line = clean_line.replace('A-', 'A -')
            clean_line = clean_line.replace('C-', 'C -')
               #feedrate is always missing a whitespace character, so add one
            clean_line = clean_line.replace('F', 'F ')
               #split the list
            split_list = clean_line.split(' ')
            
            # Parsing can be confusing, so run line-by-line and use flags to help define action for that based on what's come before
               #Ignore lines that are empty; avoids errors
            if len(line)>0:
                
                # Pull useful info from
                if toolpath_flag and (not footer_flag) and ('FREERUN' not in line):
                    line_tuple = [None, None, None, None, None, None, None]
                    update_idx = None
                    # Run through each part of the line 
                    for string in split_list:
                        #if 'string' is X,Y,Z,A,C, or F, 
                        try:
                            update_idx = idx_dict[string]
                            #if the string is 'G1' or some other mode, add to mode
                            if update_idx==6 and line_tuple[6]==None:
                                line_tuple[update_idx] = string
                            elif type(line_tuple[6])=='list':
                                line_tuple[6] = line_tuple[6].append(string)
                            elif line_tuple[6]==None:
                                line_tuple[6] = [string]
                        except:
                            pass
        
                        #if this 'string' can be converted into a float, assume the last flag is where that number should go
                        try:
                            number = float(string)
                            line_tuple[update_idx] = number
                        except:
                            pass
                    xyzacf_mode_array.append(tuple(line_tuple))
                
                # 'not toolpath_flag' meaning header lines still possible; pull header info
                elif header_flag:
                    if ';' in split_list[0]:
                        pass
                    if 'ENABLE' in line:
                        for thing in split_list[1::]:
                            pass
                    if 'Offsets' in line:
                        # try to parse each coordinate from the 'offset' line
                        #     should continue working unless prompts change
                        try:
                            part.offsets.update({split_list[1]:split_list[3]})
                        except:
                            pass
                        try:
                            part.offsets.update({split_list[4]:split_list[6]})
                        except:
                            pass
                        try:
                            part.offsets.update({split_list[7]:split_list[9]})
                        except:
                            pass
                    if 'Registration' in line:
                        # Will fail last few for 3-axis/flat plates
                        #    Susceptible if ANY changes to registration prompts occurs
                        try:
                            part.registration_points.update({split_list[2]:split_list[3]})
                        except:
                            pass
                        try:
                            part.registration_points.update({split_list[4]:split_list[5]})
                        except:
                            pass
                        try:
                            part.registration_points.update({split_list[6]:split_list[7]})
                        except:
                            pass
                        try:
                            part.registration_points.update({split_list[8]:split_list[9]})
                        except:
                            pass
                        try:
                            part.registration_points.update({split_list[10]:split_list[11]})
                        except:
                            pass
                    
                # Raise flags for parsing starting next line
                if '-- HEADER --' in line:
                    header_flag = True
        
                if 'MAIN TOOL PATH' in line:
                    toolpath_flag = True
        
                if '-- FOOTER --' in line:
                    footer_flag = True
        
        return xyzacf_mode_array


    #@staticmethod
    def show_config_tree(self):
        dict = self.config_dict
        self.recursive_dictkey_printer(dict, 0)
