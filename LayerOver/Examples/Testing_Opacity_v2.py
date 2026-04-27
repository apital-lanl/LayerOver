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

#Size of the resulting generation array
array_dims = (1000,1000)
#Number of "identical" structures to generate
number_of_duplicates = 7
#Structure to generate
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
        'mech_data_flag': None,
        'keyence_data_flag': None,
        'punch_diameter': None,
        'mass_g': None,
        'thickness_mm': None,
        'density_g/cc': None
        },
    'part_structure': 'S8HS',
    'part_skin_nozzle_size': 150,
    'part_layer_nozzle_size': 150,
    'part_angular_offset': 45,
    'part_lateral_offset': 0,
    'part_material': '',
    'part_pitch': None,
    'number_of_layers': 10,
    'layer_strand_extrusion': None,
    'layer_strand_diameter': [150,150,150,150,150,150,150,150,150,150],
    'layer_types': ['skin', 'helicoidal', 'helicoidal', 'helicoidal', 'helicoidal', 'helicoidal', 'helicoidal', 'helicoidal', 'helicoidal', 'skin'],
    'layer_type_modifiers': ['None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None'],
    'layer_points': [[],[],[],[],[],[],[],[],[],[]],
    'layer_steps': [0, 125, 250, 250, 250,250,250,250,250,250],
    'layer_angles': [0, 45, 90, 135, 180, 225, 270, 315, 0, 45],
    'layer_lateral_offsets': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    'layer_materials': ['LL50', 'LL50', 'LL50', 'LL50', 'LL50', 'LL50', 'LL50', 'LL50', 'LL50', 'LL50'],
    'layer_pitches': [450, 450, 450, 450, 450, 450, 450, 450, 450, 450],
    'layer_height_modifiers': [0.78 for i in range(10)],
    'layer_heights': [0, 125, 250, 250, 250,250,250,250,250,250]
    }

#Do the generatin'
for i in range(number_of_duplicates):
    print()
    print('#'*50)
    print(f"Running replicate {i+1}")
    try:
        flat_ideal_volume_guess(diw_structure_dict,
                                    array_dims, 
                                    n_structures= 1, 
                                    voxel_side_length= 15,
                                    voxel_resolution_microns= 0,
                                    compression_factor= 0.7,
                                    save_layer_images= False,
                                    save_layer_arrays= False,
                                    save_final_image= False,
                                    save_final_array= False,
                                    show_layer_images= False,
                                    show_final_image= True)
    except Exception as e:
        print()
        print(f'Failed on exception: {e}')
        print(f"\t {traceback.format_exc()}")


#Returned keys:
#     metadata
#     number_of_layers
#     layer_strand_extrusion
#     layer_strand_diameter
#     layer_types
#     layer_type_modifiers
#     layer_points
#     layer_steps
#     layer_angles
#     layer_lateral_offsets
#     layer_materials
#     layer_pitches