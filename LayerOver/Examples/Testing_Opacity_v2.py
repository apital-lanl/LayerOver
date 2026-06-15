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

import csv
import numpy as np
import math
import matplotlib.pyplot as plt
import os
from scipy.signal import find_peaks
import traceback
from tkinter import filedialog, Tk

from LayerOver.Analysis.VolumetricPrediction import flat_ideal_volume_guess
from LayerOver.PSPP.DIWStructure import blank_diw_structure_dict
from LayerOver.PSPP.DIWStructure import fill_in_structure_dict
from LayerOver.Analysis.ImageData import dynamic_threshold

###############################################################################################################
###  User parameters  #########################################################################################
###############################################################################################################

#Main user options
  # individual layer options
show_each_layer= False
show_layer_overlaps = True
show_layer_adjustment_diffs = False
  # final 'stacked' part images
show_full_prediction= True
show_full_overlap = True
  # show analysis options
show_prediction_histograms = True
  # save/store
save_histogram = False

#Other options
overwrite_layer_lists = True

#Structure parameters
  #size of the resulting generation array
array_dims = (1000,1000)
  #number of "identical" structures to generate
number_of_duplicates = 1
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
    'part_structure': 'S5H',
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

if save_histogram:
    #Create save directory if it doesn't exist
    root = Tk()
    save_dir = filedialog.askdirectory(title='Select directory to save histogram CSVs')
    root.destroy()

    #Create a part_name string
    structure = diw_structure_dict['part_structure']
    skin_strand = diw_structure_dict['part_skin_nozzle_size']
    layer_strand = diw_structure_dict['part_layer_nozzle_size']
    angle = diw_structure_dict['part_angular_offset']
    lateral = diw_structure_dict['part_lateral_offset']
    pitch = diw_structure_dict['part_pitch']
    part_name = f"{structure}_skin({skin_strand})_layer({layer_strand})_angle({angle})_lat({lateral})_pitch({pitch})"

    #Intialize a dictionary to hold histogram data for each volume prediction
    volume_histogram_dict = {}


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
        layer_overlaps = volume_dict[volume_key]['layer_overlaps']
        adjusted_layers = volume_dict[volume_key]['adjusted_layers']
        ideal_layers = volume_dict[volume_key]['ideal_layers']

        #Returned keys for each 'volume_dict[voxel_key]':
        # 'ideal_layers'            list; each entry is a list of 2D arrays, one per layer, with ideal layer predictions
        # 'adjusted_layers'         list; each entry is a list of 2D arrays, one per layer, with adjusted layer predictions
        # 'layer_overlaps'          list; thickness of 'overlap' for each layer (in microns)
        # 'full_volume_prediction'  np.array; 2D array of stacked 'adjusted_layers' arrays
        # 'full_ideal_prediction'   np.array; 2D array of stacked 'ideal_layers' arrays
        # 'full_overlap_prediction' np.array; 2D array of stacked adjustments made to layers to account for overlap (i.e. the difference between 'full_volume_prediction' and 'full_ideal_prediction')

        if show_full_overlap:
            plt.imshow(overlap_array)
            plt.title(f"Full overlap for iteration {i}")
            plt.show()
        
        if show_layer_overlaps:
            for idx, overlap_array in enumerate(layer_overlaps):
                plt.imshow(overlap_array)
                plt.title(f"Layer overlap array for iteration {i}-layer {idx+1} of {len(layer_overlaps)}")
                plt.show()

        if show_layer_adjustment_diffs:
            for idx, (this_ideal_array, adjusted_array) in enumerate(zip(ideal_layers, adjusted_layers)):
                this_diff_array = this_ideal_array-adjusted_array
                plt.imshow(this_diff_array)
                plt.title(f"Difference between ideal and compressed; iteration {i}-layer{idx+1} of {len(ideal_layers)}")
                plt.show()

        if show_prediction_histograms or save_histogram:
            comp_hist_dict = dynamic_threshold(volume_array, 
                            num_bins = 200,
                            show_threshold_graph = show_prediction_histograms,
                            name = f'Duplicate {i} compression factored prediction', 
                            threshold = True,
                            calculation_range = 'below',
                            calculation_type = 'outside_CLT',
                            fix_bins = True)

                # ideal_hist_dict keys:
                #     'name'                    str
                #     'cnts'                    list; histogram counts
                #     'bins'                    list; bin edges
                #     'max_middle_cnt'          int; maximum count in the middle bins
                #     'bin_size'                float; size of each bin
                #     'lower_FWHM_bin_idx'      int; index of the lower half-maximum bin
                #     'upper_FWHM_bin_idx'      int; index of the upper half-maximum bin
                #     'FWHM'                    float; full width at half maximum counts
                #     'calculation'             str; calculation type and range
                #     'calculation_cnt_sum'     int; sum of counts used in calculation
                #     'upper_FWHM_pix_value'    float; pixel value at upper half-maximum
                #     'peak_pix_value'          float; pixel value at peak
                #     'lower_FWHM_pix_value'    float; pixel value at lower half-maximum
                #     'max_threshold_pix_value' float; maximum threshold pixel value
                #     'min_threshold_pix_value' float; minimum threshold pixel value

            # dynamic_threshold(ideal_array, 
            #                 num_bins = 200,
            #                 show_threshold_graph = True,
            #                 name = f'Duplicate {i} ideal prediction', 
            #                 threshold = True,
            #                 calculation_range = 'below',
            #                 calculation_type = 'outside_CLT',
            #                 fix_bins = True)

            overlap_hist_dict = dynamic_threshold(overlap_array, 
                            num_bins = 100,
                            show_threshold_graph = True,
                            name = f'Duplicate {i} compression compensation', 
                            threshold = True,
                            calculation_range = 'below',
                            calculation_type = 'outside_CLT',
                            fix_bins = True)

            comp_peak_bins = comp_hist_dict['bins'][1::]
            comp_peak_cnts = comp_hist_dict['cnts']
            overlap_bins = overlap_hist_dict['bins'][1::]
            overlap_cnts = overlap_hist_dict['cnts']
            comp_peak_indcs, comp_peak_props = find_peaks(comp_peak_cnts)

            #Create combined spatial overlap-thickness histogram
              # flatten each 2D grid to a 1D sample vector, then make N×2 array
            stacked_array = np.column_stack((volume_array.ravel(), overlap_array.ravel())) 
              # compute the 3D joint histogram
            counts_2d, edges_2d = np.histogramdd(stacked_array, bins=(100, 100))
            #'edges_2d'- [0] = thickness edges
            #            [1] = 
              # plot the results
            # plt.imshow(counts_2d, extent=[edges_2d[0][0],edges_2d[0][-1],edges_2d[1][0],edges_2d[1][-1]])
            plt.imshow(counts_2d)
            plt.title(f"Co-Location Histogram of Overlap and Thickness")
            plt.xlabel("Opacity values")
            plt.ylabel("Thickness values")
            plt.show()
            
            if show_prediction_histograms:
                plt.plot(comp_peak_bins[:-1], comp_peak_cnts)
                for peak_idx in comp_peak_indcs:
                    plt.scatter(comp_peak_bins[peak_idx+1], comp_peak_cnts[peak_idx], color='red', marker='x')
                plt.show()

            if save_histogram:
                volume_histogram_dict.update({f'thickness-bins, replicate{i}': comp_peak_bins})
                volume_histogram_dict.update({f'thickness-cnts, replicate{i}': comp_peak_cnts})
                volume_histogram_dict.update({f'overlap-bins, replicate{i}': overlap_bins})
                volume_histogram_dict.update({f'overlap-cnts, replicate{i}': overlap_cnts})
               
    except Exception as e:
        print()
        print(f'Failed on exception: {e}')
        print(f"\t {traceback.format_exc()}")


if save_histogram:
    #Save 'volume_histogram_dict' as a CSV no matter what
    savedict_filepath = os.path.join(save_dir, f'volume_histogram_{part_name}.csv')
    with open(savedict_filepath, 'w', newline='') as file:
        writer = csv.writer(file)
        for key, value in volume_histogram_dict.items():
            writer.writerow([key, value])


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