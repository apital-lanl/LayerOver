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
#Standard libraries
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import skimage.draw as draw

#LayerOver imports
  #separate calls for clarity; 
from LayerOver.Core.Points import draw_2D_strand_line_bythickness
from LayerOver.Core.Points import generate_random_pole_point_2D
from LayerOver.Core.Points import pole_point_to_array_interior_point
from LayerOver.Core.Points import ideal_tile_from_initial_line

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
###  Volume production from structure data  #####################################################################
#################################################################################################################

def flat_ideal_volume_guess(structure_dict,
                            array_dims, 
                            n_structures= 1, 
                            voxel_side_length= 0,
                            voxel_resolution_microns= 0,
                            compression_factor= 0.7,
                            save_layer_images= False,
                            save_layer_arrays= False,
                            save_final_image= False,
                            save_final_array= False,
                            save_location = None,
                            show_layer_images= True,
                            show_final_image= True):
    '''
    Description: Take 'ideal' strucures 

    INPUT:
        'structure_dict'            lorem
      (optional)
        'n_structures'              Number of voxels to generate
        'voxel_side_length'         In mm; total X-Y side length of the full voxel
        'voxel_resolution_microns'  In microns (um); pixel side length 
        'compression_factor'        Percentage of 'ideal' thickness 'actual' part will be do to strand-strand convergence pre-curing
    ACTIONS:
    OUTPUTS:
        'return_volume_dict'-       Dictionary with same structure (keys) as 'generic_volume_dict'; represents
                                    results of simulating layer feature 
    '''

    #Initialize variables

      # hard-coded parameters
    figure_width = 5   #inches; for plotting and saving images
    figure_height = 5   #inches; for plotting and saving images
    figure_dpi = 600   #dpi for figure save

      # refine length scale variables
    if voxel_side_length == 0 and voxel_resolution_microns == 0:
        voxel_side_length = 15.875  #mm; default to punch diameter of 15.875 mm (5/8 inch) if no other length options are specified
    elif voxel_side_length == 0:

        #'array_dims' must be specified, so if 'voxel_resolution_microns' is specified, use it to constrain mm size scale
        voxel_side_length = max(array_dims) * voxel_resolution_microns /1000  #pixels->microns->mm
    print(f"Array side length: {voxel_side_length} mm")
    if voxel_resolution_microns == 0:
        voxel_resolution_microns = round((voxel_side_length*1000)/max(array_dims), 5)
    print(f"Pixel resolution: {voxel_resolution_microns} microns")
    um_to_pix_conversion = voxel_resolution_microns   #Renaming for clarity and ease of use in functions below
    x_array_dim= int(round(voxel_side_length *1000 / voxel_resolution_microns))   #number of pixels (microns/microns)
    y_array_dim= x_array_dim  #Array is square
    array_dims = (y_array_dim, x_array_dim)
    print(f"Array dimensions: {array_dims}")
    return_volume_dict = {
        'array_dims': array_dims}
    pixel_half_length = voxel_resolution_microns/2
    x_distance_idcs = np.linspace(pixel_half_length, (voxel_side_length *1000)-pixel_half_length, x_array_dim-1)
    y_distance_idcs = np.linspace(pixel_half_length, (voxel_side_length *1000)-pixel_half_length, y_array_dim-1)
    # mesh_x, mesh_y = np.array(np.meshgrid(x_distance_idcs, y_distance_idcs))

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

    #Check initial data types and coerce
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
        print(f"Generating structure {i}")
        #Generate a layer name
        voxel_idx = i+1
        if metadata_dict['unique_structure_name']:
            voxel_name = metadata_dict['unique_structure_name']+f'_Voxel-{voxel_idx}'
        elif 'print_name' in metadata_dict and metadata_dict['print_name']:
            voxel_name = metadata_dict['print_name']+f'_Voxel-{voxel_idx}'
          # generate generic voxel_name if everything else fails
        else:
            #TODO: lookup unique names and structure IDs to make this more meaningful
            voxel_name = f'UNK-Structure-ID_Voxel-{voxel_idx}'
        print(f"\t Voxel name: {voxel_name}")

        voxel_volume_array = np.zeros((y_array_dim, x_array_dim))  #each dim should be the same, but keeping them separable for now in case that's not true in the future
        ideal_layer_arrays = []    #store raw strand volume array for each layer
        adjusted_layer_arrays = []    #store volume array adjusted for compression and strand-strand interaction
        #Generate each layer's thickness projection
        for layer_idx in range(n_layers):
            print(f"\t\t adding layer {layer_idx +1}")
            layer_name = f"{voxel_name}_Layer-{layer_idx+1}"
            layer_dict = {
                'ideal_layers':{},
                'adjusted_layers':{},
                'full_volume_prediction':[]}
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

            #Create layer blank
              # (Y,X) format to align with image libraries (i.e. CV2, Matplotlib, etc.)
            this_layer = np.zeros((y_array_dim, x_array_dim))

            #Populate array with lines 
            if these_layer_point_coordinates:
                #TODO: Fix this case to handle actual paths 
                #If explicit coordinates are passed, use those to draw strands
                array_dict = draw_2D_strand_line_bythickness(these_layer_point_coordinates,
                                                            this_angular_offset,
                                                            strand_diameter,
                                                            (y_array_dim, x_array_dim),
                                                            um_to_pix_conversion = um_to_pix_conversion,
                                                            length = None,
                                                            line_type = 'simple',
                                                            thickness_fcn = 'cylinder',
                                                            show_points = False,
                                                            show_array_iterations = show_layer_images,
                                                            show_final_array = False)

                drawn_array = array_dict['drawn_array']
                drawn_mask = drawn_array[drawn_array > 0]
                #If 'None' is not passed, assume the array is good and add new values to 'this_layer'
                if drawn_array:
                    this_layer[drawn_mask] = drawn_array[drawn_mask]

            #If no explicit coordinates are passed, assume this is a generic layer and populate with strands as appropriate
            else:
                #For first layer, generate seed points
                if layer_idx == 0:
                    #Create a seed point 
                    seed_y_idx = random.randrange(0, y_array_dim-1)
                    seed_x_idx = random.randrange(0, x_array_dim-1)
                    starting_interior_point = [seed_y_idx, seed_x_idx]

                    #Draw initial line
                    array_dict = draw_2D_strand_line_bythickness(starting_interior_point,
                                                                this_angular_offset,
                                                                strand_diameter,
                                                                (y_array_dim, x_array_dim),
                                                                um_to_pix_conversion = um_to_pix_conversion,
                                                                length = None,
                                                                line_type = 'simple',
                                                                thickness_fcn = 'cylinder',
                                                                show_points = False,
                                                                show_array_iterations = False,
                                                                show_final_array = False)

                    #Extend initial seed point to an arbitrary 'pole_point' that's used for every other layer
                    #NOTE: 'pole_point' is the point at which the strand will be drawn through for every layer.
                    #       It is generally outside the drawn array. 
                    pole_point = generate_random_pole_point_2D([seed_y_idx, seed_x_idx], this_angular_offset, array_dims)
                
                #For each additional layer, use pole_point to generate a new line
                else:
                    #Get new layer points from 'pole_point' and new layer-structure specification
                    starting_interior_point = pole_point_to_array_interior_point(pole_point,
                                                                                this_angular_offset,
                                                                                this_pitch,
                                                                                this_lateral_offset,
                                                                                array_dims,
                                                                                strand_diameter,
                                                                                um_to_pix_conversion = um_to_pix_conversion,
                                                                                line_type = 'simple',
                                                                                )

                    #Draw initial line
                    array_dict = draw_2D_strand_line_bythickness([starting_interior_point[0], starting_interior_point[1]],
                                                                this_angular_offset,
                                                                strand_diameter,
                                                                (y_array_dim, x_array_dim),
                                                                um_to_pix_conversion = um_to_pix_conversion,
                                                                length = None,
                                                                line_type = 'simple',
                                                                thickness_fcn = 'cylinder',
                                                                show_points = False,
                                                                show_array_iterations = False,
                                                                show_final_array = False)

                #Pull array and values from this layer's run
                line_key = list(array_dict['line_dicts'].keys())[0]
                line_dict = array_dict['line_dicts'][line_key]  #Should only be one entry, so taking first value is good enough
                initial_line_start = line_dict['start_coordinates']
                initial_line_end = line_dict['end_coordinates']
                drawn_array = line_dict['drawn_array']
                drawn_mask = drawn_array > 0
                
                #If 'drawn_array' has a value that isn't None, assume the array is good and add new values to 'this_layer'
                if drawn_array is None:
                    #TODO: add error handling for this case
                    pass  
                else:
                    this_layer[drawn_mask] = drawn_array[drawn_mask]
                    

                #Add all the other lines
                tile_dict = ideal_tile_from_initial_line([initial_line_start, initial_line_end],
                                                        this_angular_offset, 
                                                        this_pitch,
                                                        array_dims,
                                                        strand_diameter,
                                                        um_to_pix_conversion= um_to_pix_conversion,
                                                        show_final_array = False)
                tile_array = tile_dict['drawn_array']
                tile_mask = tile_array > 0
                this_layer[tile_mask] = tile_array[tile_mask]
            
                #Store a 'raw' version of the layer
                ideal_layer_arrays.append(this_layer)

                #Adjust for compression, strand-to-strand interactions and add to global volume array
                if layer_idx == 0:
                    #Apply flat-plate compression (i.e. compression of strand against plate surface)
                    pass
                else:
                    last_layer_array = ideal_layer_arrays[layer_idx-1]
                    
            
                #Handle images for each layer (ideal, no compression)
                if save_location:
                    layer_save_name = os.path.join(save_location, layer_name)
                else:
                    layer_save_name = layer_name
                if show_layer_images:
                    plt.figure(figsize=(figure_width,figure_height))
                    plt.imshow(this_layer)
                    plt.title(layer_name)
                    plt.show()

                if save_layer_images:
                    fig = plt.figure(frameon=False)
                    fig.set_size_inches(figure_width,figure_height)
                    ax = plt.Axes(fig, [0., 0., 1., 1.])
                    ax.set_axis_off()
                    fig.add_axes(ax)
                    ax.imshow(this_layer, aspect='auto')
                    fig.savefig(layer_save_name, dpi = figure_dpi)

                #Save array if applicable
                if save_layer_arrays:
                    if save_location:
                        layer_array_savename = os.path.join(save_location, f"{layer_name}_RawArray")
                    else:
                        layer_array_savename = f"{layer_name}_RawArray"
                    np.save(layer_array_savename, this_layer)

                #Store the layer arrays in memory
                adjusted_layer_arrays.append(this_layer)    #save adjusted layer; accounts for compression and strand-strand interactions
                voxel_volume_array = voxel_volume_array + this_layer    #add adjusted layer to global volume array
        
            #Assign to dict for return
            layer_dict['ideal_layers'][layer_idx] = ideal_layer_arrays
            layer_dict['adjusted_layers'][layer_idx] = adjusted_layer_arrays
            layer_dict['full_volume_prediction'] = voxel_volume_array
            return_volume_dict[voxel_name] = layer_dict

        #Save and/or show results
        if save_location:
            image_name = f"FullStack_Image_{voxel_name}"
            final_image_savename = os.path.join(save_location, image_name)
        else:
            final_image_savename = f"FullStack_Image_{voxel_name}"
        if show_final_image:
            plt.figure(figsize=(10,10))
            plt.imshow(voxel_volume_array)
            plt.title(image_name)
            plt.show()
        if save_final_image:
            fig = plt.figure(frameon=False)
            fig.set_size_inches(figure_width,figure_height)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(voxel_volume_array, aspect='auto')
            fig.savefig(final_image_savename, dpi = figure_dpi)

            # save the array as a numpy file
        if save_final_array:
            if save_location:
                final_array_savename = os.path.join(save_location, f"FullStack_Array_{voxel_name}")
            else:
                final_array_savename = f"FullStack_Array_{voxel_name}"
            np.save(final_array_savename, voxel_volume_array)
        
    return return_volume_dict


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


#################################################################################################################
###  Helper utilities  ##########################################################################################
#################################################################################################################


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


#################################################################################################################
###  Analysis and Statistics  ###################################################################################
#################################################################################################################


def get_volumetric_array_statistics(array,
                                 n_histogram_bins = 50,
                                 n_rounding_places= 4):
    '''
    Description: 
        Standardized reporting function for an array of volumetric data.
    
    INPUT:
        ''      list or numpy.array of volumetric predictions

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