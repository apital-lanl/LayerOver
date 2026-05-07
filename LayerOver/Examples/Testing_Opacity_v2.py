"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-03-31
Version:   1.0.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Quick script to run 'VolumetricPrediction' generation of layer images from ideal structures.

"""

import numpy as np
import math
import matplotlib.pyplot as plt
import traceback

from LayerOver.Analysis.VolumetricPrediction import flat_ideal_volume_guess
from LayerOver.PSPP.DIWStructure import blank_diw_structure_dict
from LayerOver.PSPP.DIWStructure import fill_in_structure_dict
from LayerOver.Analysis.ImageData import dynamic_threshold

###############################################################################################################
###  User parameters  #########################################################################################
###############################################################################################################

#Main user options
show_prediction_histograms = False
show_each_layer= True
show_full_prediction= False
overwrite_layer_lists = False

#Structure parameters
  #size of the resulting generation array
array_dims = (1000,1000)
  #number of "identical" structures to generate
number_of_duplicates = 3
  #structure to generate
diw_structure_dict = {
    'metadata': {
        'unique_structure_name': None,
        'print_name': None,
        'structure': None,
        'nozzle_size_um': None,
        'pitch_offset': None,
        'syringe-material': None,
        'project': None,
        'machine_name': None,
        'layerup_file': None,
        'version_number': None,
        'notes': None,
        'mech_data_note': None,
        'mech_data_filepaths': None,
        'keyence_data_note': None,
        'punch_diameter': None,
        'mass_g': None,
        'thickness_mm': None,
        'density_g/cc': None
        },
    'part_structure': 'S8HS',
    'part_skin_nozzle_size': 150,
    'part_layer_nozzle_size': 150,
    'part_angular_offset': 40,
    'part_lateral_offset': 0,
    'part_material': '',
    'part_pitch': 500,
    'number_of_layers': None,
    'layer_strand_extrusion': None,
    'layer_strand_diameter': None,
    'layer_types': None,
    'layer_type_modifiers': None,
    'layer_points': None,
    'layer_steps': None,
    'layer_angles': None,
    'layer_lateral_offsets': None,
    'layer_materials': None,
    'layer_pitches': None,
    'layer_height_modifiers': None,
    'layer_heights': None
    }

###############################################################################################################
###  Prediction and Visualization  ############################################################################
###############################################################################################################

#Fill in structure dict with lists for each layer
diw_structure_dict= fill_in_structure_dict(diw_structure_dict,
                                           overwrite_layer_lists = True)

#Do the generatin'
for i in range(number_of_duplicates):
    print()
    print('#'*50)
    print(f"Running replicate {i+1}")
    try:
        volume_dict = flat_ideal_volume_guess(diw_structure_dict,
                                    array_dims, 
                                    n_structures= 1, 
                                    voxel_side_length= 15,
                                    voxel_resolution_microns= 0,
                                    compression_factor= 0.7,
                                    save_layer_images= False,
                                    save_layer_arrays= False,
                                    save_final_image= False,
                                    save_final_array= False,
                                    show_layer_images= show_each_layer,
                                    show_final_image= show_full_prediction)

        volume_key = list(volume_dict.keys())[0]
        volume_array = volume_dict[volume_key]['full_volume_prediction']
        ideal_array = volume_dict[volume_key]['full_ideal_prediction']
        overlap_array = volume_dict[volume_key]['full_overlap_prediction']
        
        #Returned keys for each 'layer_dict':
        # 'ideal_layers'            list; each entry is a list of 2D arrays, one per layer, with ideal layer predictions
        # 'adjusted_layers'         list; each entry is a list of 2D arrays, one per layer, with adjusted layer predictions
        # 'full_volume_prediction'  np.array; 3D array of stacked 'adjusted_layers' arrays
        # 'full_ideal_prediction'   np.array; 3D array of stacked 'ideal_layers' arrays
        # 'full_overlap_prediction' np.array; 3D array of stacked adjustments made to layers to account for overlap (i.e. the difference between 'full_volume_prediction' and 'full_ideal_prediction')

        if show_prediction_histograms:
            dynamic_threshold(volume_array, 
                            num_bins = 200,
                            show_threshold_graph = True,
                            name = f'Duplicate {i} compression factored prediction', 
                            threshold = True,
                            calculation_range = 'below',
                            calculation_type = 'outside_CLT',
                            fix_bins = True)

            dynamic_threshold(ideal_array, 
                            num_bins = 200,
                            show_threshold_graph = True,
                            name = f'Duplicate {i} ideal prediction', 
                            threshold = True,
                            calculation_range = 'below',
                            calculation_type = 'outside_CLT',
                            fix_bins = True)

            dynamic_threshold(overlap_array, 
                            num_bins = 100,
                            show_threshold_graph = True,
                            name = f'Duplicate {i} compression compensation', 
                            threshold = True,
                            calculation_range = 'below',
                            calculation_type = 'outside_CLT',
                            fix_bins = True)

    except Exception as e:
        print()
        print(f'Failed on exception: {e}')
        print(f"\t {traceback.format_exc()}")



# blank_diw_structure_dict = {
#     'metadata': {
#         'unique_structure_name': None,
#         'print_name': None,
#         'structure': None,
#         'nozzle_size_um': None,
#         'pitch_offset': None,
#         'syringe-material': None,
#         'project': None,
#         'machine_name': None,
#         'layerup_file': None,
#         'version_number': None,
#         'notes': None,
#         'mech_data_note': None,
#         'mech_data_filepaths': None,
#         'keyence_data_note': None,
#         'punch_diameter': None,
#         'mass_g': None,
#         'thickness_mm': None,
#         'density_g/cc': None
#         },
#     'part_structure': None,
#     'part_skin_nozzle_size': None,
#     'part_layer_nozzle_size': None,
#     'part_angular_offset': None,
#     'part_lateral_offset': None,
#     'part_material': '',
#     'part_pitch': None,
#     'number_of_layers': None,
#     'layer_strand_extrusion': None,
#     'layer_strand_diameter': None,
#     'layer_types': None,
#     'layer_type_modifiers': None,
#     'layer_points': None,
#     'layer_steps': None,
#     'layer_angles': None,
#     'layer_lateral_offsets': None,
#     'layer_materials': None,
#     'layer_pitches': None,
#     'layer_height_modifiers': None,
#     'layer_heights': None
#     }