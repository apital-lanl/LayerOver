"""
Copyright 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2025-11-25
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Module for handling conversions between ideal structure assumptions and physical measurements of as-printed structures.
    NOTE: Code below is for specific nomenclature and testing conventions. An attempt has been made to highlight hard-coded assumptions 
          and segregate these to the headers of classes and the full module. Other (specific) use cases will require modifications.

"""
#from Points import
import numpy as np
import random
import skimage.draw as draw
from LayerOver.Core.Points import draw_2D_strand_line_bythickness


#Hard-coded assumptions and conversion ratios

default_compression_factors= {
    'default': 0.7,
    'general-siloxane': 0.7,
    }
#Generic dictionaries 
generic_punch_sizes = {
    'default':{
        'diameter_units': 'microns',
        'diameter': 15.875
        },
    'generic_mech_testing':{
        'diameter_units': 'inches',
        'diameter': 5/8
        }
    }

# 'print_type'-  'flat'; 'hemi'; 'shaped-special'
generic_volume_dict = {
    'print_type': '',
    'number_of_layers': 0,
    'layer_types':[],
    'strand_sizes': [],
    'layer_arrays': []
    }


#################################################################################################################
###  Lorem  #####################################################################################################
#################################################################################################################

def flat_ideal_volume_guess(structure_dict, 
                       n_structures= 2, 
                       voxel_side_length= 15.875,
                       voxel_resolution_microns= 1,
                       compression_factor= 0.7,
                       save_layer_images= False,
                       save_layer_arrays= False,
                       save_final_image= False,
                       save_final_array= False,
                       show_layer_images= True,
                       show_final_image= True):
    '''
    Description: Take 

    INPUT:
        'structure_dict'-           lorem
      (optional)
        'n_structures'-             Number of voxels to generate
        'voxel_side_length'-        in mm; total X-Y side length of the full voxel
        'voxel_resolution'-         In microns (um); pixel side length 
        'compression_factor'-       Percentage of 'ideal' thickness 'actual' part will be do to strand-strand convergence pre-curing
    ACTIONS:
    OUTPUTS:
        'return_volume_dict'-       Dictionary with same structure (keys) as 'generic_volume_dict'; represents
                                    results of simulating layer feature 
    '''

    #Initialize variables
    return_volume_dict = {}
    x_array_dim= int(round(voxel_side_length *1000 / voxel_resolution_microns))   #number of pixels (microns/microns)
    y_array_dim= x_array_dim  #Array is square
    pixel_half_length = voxel_resolution_microns/2
    x_distance_idcs = np.linspace(pixel_half_length, (voxel_side_length *1000)-pixel_half_length, x_array_dim-1)
    y_distance_idcs = np.linspace(pixel_half_length, (voxel_side_length *1000)-pixel_half_length, y_array_dim-1)
    mesh_x, mesh_y = np.array(np.meshgrid(x_distance_idcs, y_distance_idcs))
      # pull 'structure_dict' keys
    metadata_dict= structure_dict['metadata']
    n_layers= structure_dict['number_of_layers']
    layer_extrusion_type= structure_dict['layer_strand_extrusion']
    strand_diameters= structure_dict['layer_strand_diameter']
    layer_types= structure_dict['layer_types']
    layer_modifiers= structure_dict['layer_type_modifiers']
    layer_point_coordinates= structure_dict['layer_points']
    layer_heights= structure_dict['layer_steps']
    layer_angular_offsets= structure_dict['layer_angles']
    layer_lateral_offsets= structure_dict['layer_lateral_offsets']
    layer_materials= structure_dict['layer_materials']
    layer_pitches= structure_dict['layer_pitches']

    #Check data types and coerce
    if type(compression_factor) == float:
        # check if a fraction
        if compression_factor <=1 and compression_factor >0:
            pass
        else:
            #TODO: add flag for non-fractional compression factors
            pass
    elif type(compression_factor) == int:
        if compression_factor >1 and compression_factor <=100:
            compression_factor = round(compression_factor/100, 5)
        elif compression_factor == 1:
            #Assume no compression and just pass as 1; otherwise compression will be unphysical
            pass
    
    #Generate 'n_structures' number of iterations of volume
    for i in range(n_structures):
        
        #Generate a layer name
        voxel_idx = i+1
        if 'unique_structure_name' in metadata_dict:
            voxel_name = metadata_dict['unique_structure_name']+f'_Voxel-{voxel_idx}'
        elif 'print_name' in metadata_dict:
            voxel_name = metadata_dict['print_name']+f'_Voxel-{voxel_idx}'
          # generate generic voxel_name if everything else fails
        else:
            #TODO: lookup unique names and structure IDs to make this more meaningful
            voxel_name = f'UNK-Structure-ID_Voxel-{voxel_idx}'

        voxel_volume_array = np.zeros((y_array_dim, x_array_dim))  #each dim should be the same, but keeping them separable for now in case that's not true in the future
        layer_arrays = []
        #Generate each layer's thickness projection
        for layer_idx in range(n_layers):
            #Get layer specifics
            strand_diameter = strand_diameters[layer_idx]
            this_layer_type = layer_types[layer_idx]
            this_layer_modifier = layer_modifiers[layer_idx]
            these_layer_point_coordinates = layer_point_coordinates[layer_idx]
            this_layer_height = layer_heights[layer_idx]
            this_angular_offset = layer_angular_offsets[layer_idx]
            this_lateral_offset = layer_lateral_offsets[layer_idx]
            this_material = layer_materials[layer_idx]
            this_pitch = layer_pitches[layer_idx]
              # convert strand_diameter to # of pixels
            strand_radius_in_pixels = round(strand_diameter/2/voxel_resolution_microns, 5)
            

            #Get layer matrial opacity 

            #Create layer blank
              # (Y,X) format to align with image libraries (i.e. CV2, Matplotlib, etc.)
            this_layer = np.zeros((y_array_dim, x_array_dim))

            #Populate array with lines 
            if these_layer_point_coordinates:
                #If explicit coordinates are passed, use those to draw strands
                draw_2D_strand_line_bythickness(these_layer_point_coordinates,
                                                this_angular_offset,
                                                strand_radius_in_pixels,
                                                (y_array_dim, x_array_dim),
                                                length = None,
                                                line_type = 'simple')
            #If no explicit coordinates are passed, assume this is a generic layer and populate with strands as appropriate
            else:
                #Create a seed point 
                seed_y_idx = random.randrange(0, y_array_dim-1)
                seed_x_idx = random.randrange(0, x_array_dim-1)

            #Populate array with strand thicknesses
              # find edge coordinates for the seed point


              # draw initial line based on seed

              # check max distances and fill in the rest of the array
            
            #Adjust for compression, strand-to-strand interactions and add to global volume array
            if layer_idx == 0:
                pass
            else:
                pass

            if save_layer_arrays:
                pass
            
            #Handle images for each layer (ideal, no compression)
            if show_layer_images:
                if save_layer_images:
                    pass
                
              # save the layer image if 'show_layer_images' is False
            elif save_layer_images:
                pass

            layer_arrays.append(this_layer)

        #



def cartesian_ideal_volume_guess(structure_dict, 
                       n_structures= 5, 
                       voxel_side_length= 15.875,
                       voxel_resolution_microns= 1,
                       compression_factor= 0.7,
                       save_layer_images= False,
                       save_layer_arrays= False,
                       save_final_image= False,
                       save_final_array= False,
                       show_layer_images= True,
                       show_final_image= True):
    '''
    Description: Take geometry from a pointcloud and generate and idealized voxel.

    INPUT:
        'structure_dict'-           lorem
      (optional)
        'n_structures'-             Number of voxels to generate
        'voxel_side_length'-        in mm; total X-Y side length of the full voxel
        'voxel_resolution'-         In microns (um); pixel side length 
        'compression_factor'-       Percentage of 'ideal' thickness 'actual' part will be do to strand-strand convergence pre-curing
    ACTIONS:
    OUTPUTS:
        'return_volume_dict'-       Dictionary with same structure (keys) as 'generic_volume_dict'; represents
                                    results of simulating layer feature 
    '''

    pass


def get_array_edges_for_line(seed_indices, 
                             array_dimensions, 
                             line_angle):
    '''
    Descscription: Get the array edge points for an arbitrary line
    '''


    pass


def apply_layer_slump(layers_array, 
                      compression_factor= None, 
                      structure_dict = None):
    """
    Description:
        Lorem
    INPUT:
        ''              lorem
    ACTION:
        -lorem
    OUTPUT:
        ''

    """

    #Initialize variables


    pass