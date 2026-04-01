"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-03-30
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: 

"""

#Import libraries and functions
import numpy as np
import os
import pandas as pd
from tkinter import Tk, filedialog
  # other LayerOver imports
from LayerOver.Analysis.MechanicalData import open_logbook
from LayerOver.Analysis.VolumetricPrediction import flat_ideal_volume_guess
from LayerOver.PSPP.DIWStructure import structure_dict_from_logbook_row

#Volume array settings
array_dims  = (5000, 5000)
number_of_replicates = 3    #Number of separate volumetric predictions to run
array_side_length_mm = 15   #size of the array in mm (everything gets scaled)
save_to_filepath = ''       #Optional; where to save images and arrays to;
    # ''                    save to same filepath that logbook is in
    # 'actual/filepath'     manual entry of valid filepath
    # 'dialog'

#Select the file to open
root = Tk()
logbook_filepath = filedialog.askopenfilename(filetypes = [('Logbook CSV Files', '*AutomatedAnalysis.csv')])
root.destroy()
  # condition a savepath to the same location if none is chosen
if save_to_filepath == '':
    save_to_filepath = os.path.dirname(logbook_filepath)
elif save_to_filepath == 'dialog':
    root = Tk()
    save_to_filepath = filedialog.askdirectory(title="Select a directory location to save images and/or arrays to.")
    root.destroy()
else:
    if not os.path.exists(save_to_filepath):
        root = Tk()
        save_to_filepath = filedialog.askdirectory(title="Select a directory location to save images and/or arrays to.")
        root.destroy()

#Open the file as a pd.DataFrame
logbook_df = open_logbook(logbook_filepath,
                          target_excel_sheetname = None)

#Iterate through data rows
print_names = logbook_df['Name']
name_list = []
structure_list = []
unique_structure_list = []
diameter_list = []
angle_list = []
lateral_list = []
pitch_list = []
thickness_list = []
density_list = []
vox_1_average = []
vox_1_sum = []
vox_1_density = []
vox_2_average = []
vox_2_sum = []
vox_2_density = []
vox_3_average = []
vox_3_sum = []
vox_3_density = []
for idx, name in enumerate(print_names):
    print('#'*50)
    print(f'Processing row {idx +1}: {name}')

    #Add name to list
    name_list.append(name)

    try:
        #Pull logbook row and generate a structure dict from it
        this_logbook_row = logbook_df[logbook_df['Name'] == name]
          # store values for later
        structure_list.append(this_logbook_row[''].values[0])
        diameter_list.append(this_logbook_row[''].values[0])
        angle_list.append(this_logbook_row[''].values[0])
        lateral_list.append(this_logbook_row[''].values[0])
        pitch_list.append(this_logbook_row[''].values[0])
        thickness_list.append(this_logbook_row[''].values[0])
        density_list.append(this_logbook_row[''].values[0])
          # generate structure dict
        diw_structure_dict = structure_dict_from_logbook_row(this_logbook_row)
    except Exception as e:
        print()
        print(f'Failed structure parsing on exception: {e}')

        structure_list.append(np.nan)
        diameter_list.append(np.nan)
        angle_list.append(np.nan)
        lateral_list.append(np.nan)
        pitch_list.append(np.nan)
        thickness_list.append(np.nan)
        density_list.append(np.nan)

    #Run the structure to try and generate a volume prediction from it
    try:
        volume_dict = flat_ideal_volume_guess(diw_structure_dict,
                                            array_dims, 
                                            n_structures= number_of_replicates, 
                                            voxel_side_length= array_side_length_mm,
                                            voxel_resolution_microns= 0,
                                            compression_factor= 0.7,
                                            save_layer_images= True,
                                            save_layer_arrays= True,
                                            save_final_image= True,
                                            save_final_array= True,
                                            save_location = save_to_filepath
                                            show_layer_images= False,
                                            show_final_image= False)
            #Returned keys:
            # <voxel_name>
            #       ideal_layers
            #       adjusted_layers
            #       full_volume_pred

        for vox_idx, voxel_key in enumerate(list(volume_dict.keys())):
            this_dict = volume_dict[voxel_key]
            this_full_volume_array = this_dict['full_volume_pred']
            #Only pulling the first three (if more are done)
            if vox_idx < 3:
                this_average = np.average(this_full_volume_array)
                this_sum = np.sum(this_full_volume_array)
                this_max = np.max(this_full_volume_array)
                max_volume = this_max * array_dims[0] * array_dims[1]
                this_fractional_density = max_volume/this_sum
                if vox_idx == 0:
                    vox_1_average.append(this_average)
                    vox_1_sum.append(this_sum)
                    vox_1_density.append(this_fractional_density)
                if vox_idx == 1:
                    vox_2_average.append(this_average)
                    vox_2_sum.append(this_sum)
                    vox_2_density.append(this_fractional_density)
                if vox_idx == 2:
                    vox_3_average.append(this_average)
                    vox_3_sum.append(this_sum)
                    vox_3_density.append(this_fractional_density)


    except Exception as e:
        print()
        print(f'Failed volume prediction on exception: {e}')

        #Add dummy variables
        vox_1_average.append(np.nan)
        vox_1_sum.append(this_sum)
        vox_1_density.append(this_fractional_density)
        vox_2_average.append(this_average)
        vox_2_sum.append(this_sum)
        vox_2_density.append(this_fractional_density)
        vox_3_average.append(this_average)
        vox_3_sum.append(this_sum)
        vox_3_density.append(this_fractional_density)




# #Save DataFrame to CSV with a similar name
# savebook_filename = os.path.basename(logbook_filepath).replace('AutomatedAnalysis', 'AutomatedResults')
# savebook_filepath = os.path.join(os.path.dirname(logbook_filepath), savebook_filename)
# with open(savebook_filepath, 'w') as file:
#     logbook_df.to_csv(file, index= False, lineterminator='\n')
