"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-03-30
Version:   1.0.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Open a logbook copy CSV (*AutomatedAnalysis.csv), create ideal structure simulations from logbook structures, and save/show/store results.

Notes:
    -Compression of strands between layers and with baseplate at first layer is accounted for in code, but is currently not implemented
    -Compressed strand layer images are saved separately in 'volume_dict'
        - 'ideal_layers' are uncompressed
        - 'adjusted_layers' have compression and other real-world accomodations made
        - 'full_volume_pred' is the entire layer stack of 'adjusted_layers'
    -If logbook entries are in rough shape, and unique structures will be created and saved for entries that make no sense (i.e. 0 radius strand layers)


"""

#Import libraries and functions
import numpy as np
import os
import pandas as pd
  # silence annoying Pandas notifications about chained assignments of row data; not an issue for this code case, but it is a bad idea generally
pd.options.mode.chained_assignment = None  # default='warn'
from tkinter import Tk, filedialog
import traceback
  # other LayerOver imports
from LayerOver.PSPP.DIWLogbook import open_logbook
from LayerOver.Analysis.VolumetricPrediction import flat_ideal_volume_guess
from LayerOver.PSPP.DIWStructure import structure_dict_from_logbook_row
from LayerOver.PSPP.DIWStructure import match_structure_to_unique_name

#Main settings
print_name_start_row  = 13   #default is 0; starting logbook row 
print_name_end_row = 0   #default is 0; ending logbook row
save_final_unique_structure_dict = False
save_to_filepath = 'W:\Data & ML\2026-04_Waffle Testing\2026-04-01_Working Alpha'       #Optional; where to save images and arrays to;
    # ''                    save to same filepath that logbook is in
    # 'actual/filepath'     manual entry of valid filepath
    # 'dialog'


#Other settings and initializations
#Volume array settings
array_dims  = (5000, 5000)
number_of_replicates = 3    #Number of separate volumetric predictions to run
array_side_length_mm = 15   #size of the array in mm (everything gets scaled)
#Initialize dict for storing structures to check for uniqueness
unique_structure_dict ={}

#Select the file to open
root = Tk()
logbook_filepath = filedialog.askopenfilename(title = "Select a logbook 'AutomatedAnalysis.csv' file", \
                                              filetypes = [('Logbook CSV Files', '*AutomatedAnalysis.csv')])
root.destroy()

#Condition a savepath to the same location if none is chosen
#NOTE: Same filepath is used to save resulting 'unique_structure_dict' if option is chosen.
#      This script assumes the layer and final volume prediction images will be saved, but that can also be turned off.
if save_final_unique_structure_dict:
      
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
if (print_name_end_row>0) and (print_name_end_row>=print_name_start_row):
    print_names = logbook_df['Name'].values[(print_name_start_row-1): (print_name_end_row)]
else:
    print_names = logbook_df['Name'].values[(print_name_start_row-1)::]

#Run through the logbook
#NOTE: If any 'print_names' are repeated, only the first version in the spreadsheet will be used
for idx, name in enumerate(print_names):
    print('#'*50)
    print(f'Processing logbook row {idx + print_name_start_row + 2}: {name}')

    #Put all of this in 'try/except' frame because logbook data is ugly and throws a ton of errors
    try:
        #Pull logbook row
        this_logbook_row = logbook_df[logbook_df['Name'] == name]
    
        #Generate structure dict
        diw_structure_dict = structure_dict_from_logbook_row(this_logbook_row)
        print(f"Structure {diw_structure_dict['part_structure']}; Strand {diw_structure_dict['part_skin_nozzle_size']}\{diw_structure_dict['part_layer_nozzle_size']}; Angle {diw_structure_dict['part_angular_offset']}; Offset {diw_structure_dict['part_lateral_offset']}")

        #Check if this is a new structure or not
        is_unique_bool, unique_structure_id, trial_unique_structure_dict = match_structure_to_unique_name(diw_structure_dict, unique_structure_dict)
        print(f"\t Structure {diw_structure_dict['part_structure']} flagged as {unique_structure_id}")

        if is_unique_bool:
            #Add some fields to this entry
            trial_unique_structure_dict[unique_structure_id].update({'example_printname': name})

            #Run the structure to try and generate a volume prediction from it
            volume_dict = flat_ideal_volume_guess(diw_structure_dict,
                                                array_dims, 
                                                n_structures= number_of_replicates, 
                                                voxel_side_length= array_side_length_mm,
                                                voxel_resolution_microns= 0,
                                                compression_factor= 0.7,
                                                save_layer_images= True,
                                                save_layer_arrays= False,
                                                save_final_image= True,
                                                save_final_array= False,
                                                save_location = save_to_filepath,
                                                show_layer_images= False,
                                                show_final_image= False)
                #Returned keys:
                # <voxel_name>
                #       ideal_layers
                #       adjusted_layers
                #       full_volume_pred

            if len(volume_dict) > 0:
                for vox_idx, voxel_key in enumerate(list(volume_dict.keys())):
                    this_dict = volume_dict[voxel_key]
                    this_full_volume_array = this_dict['full_volume_prediction']
                    #Only pulling the first three (if more are done)
                    if vox_idx < 3:
                        this_average = np.average(this_full_volume_array)
                        this_sum = np.sum(this_full_volume_array)
                        this_max = np.max(this_full_volume_array)
                        max_volume = this_max * array_dims[0] * array_dims[1]
                        this_fractional_density = this_sum/max_volume
                        if vox_idx == 0:
                            trial_unique_structure_dict[unique_structure_id].update({'vox_1_average': this_average})
                            trial_unique_structure_dict[unique_structure_id].update({'vox_1_sum': this_sum})
                            trial_unique_structure_dict[unique_structure_id].update({'vox_1_density': this_fractional_density})
                        if vox_idx == 1:
                            trial_unique_structure_dict[unique_structure_id].update({'vox_2_average': this_average})
                            trial_unique_structure_dict[unique_structure_id].update({'vox_2_sum': this_sum})
                            trial_unique_structure_dict[unique_structure_id].update({'vox_2_density': this_fractional_density})
                        if vox_idx == 2:
                            trial_unique_structure_dict[unique_structure_id].update({'vox_3_average': this_average})
                            trial_unique_structure_dict[unique_structure_id].update({'vox_3_sum': this_sum})
                            trial_unique_structure_dict[unique_structure_id].update({'vox_3_density': this_fractional_density})

                #Add specifics of how the array used in calculating voxels
                trial_unique_structure_dict[unique_structure_id].update({'full_voxel_dimensions': array_dims})
                trial_unique_structure_dict[unique_structure_id].update({'array_side_length_mm': array_side_length_mm})
                  #'array_dims' should be square, but take the max anyway; if not square, length scales will be incorrect
                trial_unique_structure_dict[unique_structure_id].update({'vox_resolution_micron': (array_side_length_mm*1000)/max(array_dims)})
                            
                #If everything worked and no errors thrown, update 'unique_structure_dict' to the new 'trial_unique_structure_dict' with the new structure
                unique_structure_dict = trial_unique_structure_dict

            #Case for when volume prediction fails in a 'normal' way (i.e. a prediction is built, but is faulty and returns None)
            else:
                print()
                print(f"Failure to build volume prediction.")
                print()
        
        # 'is_unique_bool' is False because a match was found in 'unique_structure_dict' or because the search threw errors.
        else:
            print("Not a unique structure. Moving on.")
            print()

    except Exception as e:
        print()
        print(f"Volume prediction failed with error: {e}")
        print()
        print(f"\t {traceback.format_exc()}")

        #If errors are thrown, save a rolling copy of the 'unique_structures_dict' each time as a precaution
        #NOTE: not necessary, and can be removed for speed if needed
          # format dict as DataFrame for saving
        unique_structure_df = pd.DataFrame(unique_structure_dict)
        unique_structure_df = unique_structure_df.reset_index()
          # save to location
        if save_final_unique_structure_dict:
            #Save DataFrame to CSV with a similar name
            savebook_filename = "Unique_structure_dict.csv"
            savebook_filepath = os.path.join(os.path.dirname(logbook_filepath), savebook_filename)
            with open(savebook_filepath, 'w') as file:
                unique_structure_df.to_csv(file, index= False, lineterminator='\n')

#Format dict as DataFrame for saving
unique_structure_df = pd.DataFrame(unique_structure_dict)
unique_structure_df = unique_structure_df.reset_index()
#Save to location
if save_final_unique_structure_dict:
    #Save DataFrame to CSV with a similar name
    savebook_filename = "Unique_structure_dict.csv"
    savebook_filepath = os.path.join(os.path.dirname(logbook_filepath), savebook_filename)
    with open(savebook_filepath, 'w') as file:
        unique_structure_df.to_csv(file, index= False, lineterminator='\n')
