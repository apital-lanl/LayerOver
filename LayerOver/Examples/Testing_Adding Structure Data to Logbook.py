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
pd.options.mode.chained_assignment = None  # default='warn'
from tkinter import Tk, filedialog
import traceback
  # other LayerOver imports
from LayerOver.PSPP.DIWLogbook import open_logbook
from LayerOver.PSPP.DIWStructure import open_wafflebook
from LayerOver.Analysis.VolumetricPrediction import flat_ideal_volume_guess
from LayerOver.PSPP.DIWStructure import structure_dict_from_logbook_row
from LayerOver.PSPP.DIWStructure import match_structure_to_unique_name
  #list of fields that come from PSPP analysis: mechanical data, ideal structure prediction, etc.
from LayerOver.PSPP.DIWStructure import additional_structure_fields


#######################################################################################################################
#####  Important user-defined settings  ###############################################################################
#######################################################################################################################

#Main settings
print_name_start_row  = 1   #default is 1; starting logbook row 
print_name_end_row = 0   #default is 0; ending logbook row
update_with_mechanical_data = False #option to add mech data; 
    #NOTE: will try the following in order: 1) look for scraped mech data, 2) search directory for appropriate data, 3) give up

#######################################################################################################################
#####  Main script body  ##############################################################################################
#######################################################################################################################

#Initialize variables
update_key_names = additional_structure_fields   #additional structure values not in 'diw_structure_dict' 
  # select a logbook to update
root = Tk()
logbook_filepath = filedialog.askopenfilename(title = "Select a logbook 'AutomatedAnalysis.csv' file", \
                                              filetypes = [('Logbook CSV Files', '*AutomatedAnalysis.csv')])
root.destroy()
  # select a dataset to add (from "WaffleData.csv") file
root = Tk()
waffledata_filepath = filedialog.askopenfilename(title = "Select a logbook 'WaffleData' file", \
                                              filetypes = [('Waffle Structure Summary Files', '*WaffleData.csv')])
root.destroy()
if update_with_mechanical_data:
    root = Tk()
    mechdata_directory = filedialog.askdirectory(title = "Select a directory with mechanical summary data or raw mech data")
    root.destroy()

  # open the files as a pd.DataFrames
logbook_df = open_logbook(logbook_filepath,
                          target_excel_sheetname = None)
waffledata_df = open_wafflebook(waffledata_filepath)
unique_structure_dict = waffledata_df.set_index('index').to_dict(orient='index')
unique_structureids = list(unique_structure_dict.keys())

    #'wafflebook' standard columns (i.e. columns in logbook)
        # part_structure	
        # layer_strand_diameter	
        # layer_types	
        # layer_angles	
        # layer_lateral_offsets	
        # layer_materials	
        # layer_pitches	
        # example_printname
      #Specific columns (i.e. 'wafflebook' specific from volumetric predictions)
        # index	                    str; form '<structure>_<unique structure index> (i.e. "SHS_4")
        # vox_1_average             float; average volume within the array pixels 	
        # vox_1_sum	
        # vox_1_density	
        # vox_2_average	
        # vox_2_sum	
        # vox_2_density	
        # vox_3_average	
        # vox_3_sum	
        # vox_3_density
        # full_voxel_dimensions     tuple; array dimensions used to calculated vox values above (i.e. '(1000, 1000)')
        # array_side_length_mm      int,float; length in mm of the largest array dim (should be square, so total size)
        # vox_resolution_micron     float; side-length of each pixel in the array

#Pull values from the DataFrames and condition for next steps
  # pull the unique 'print_names' to use as lookup values for each row
    # subtract 1 for 0-indexing offset 
if (print_name_end_row>0) and (print_name_end_row>=print_name_start_row):
    print_names = logbook_df['Name'].values[(print_name_start_row-1): (print_name_end_row)]
else:
    print_names = logbook_df['Name'].values[(print_name_start_row-1)::]
  # pull unique waffle structure IDs 
  # copy the logbook DataFrame to serve as the base for the updated DataFrame that's saved at the end of this
# results_df = logbook_df.copy(deep= True)
results_df = logbook_df.copy()
  # add columns to DataFrame to be updated with structure-specific data from
  #NOTE: must be in "LayerOver.PSPP.DIWStructure.additional_structure_fields" to be filled-in (updated) from 'wafflebook'
results_df['unique_structure_id'] = None   #this is the 'index' field
results_df['vox_1_average'] = None
results_df['vox_1_sum'] = None
results_df['vox_1_density'] = None
results_df['vox_2_average'] = None
results_df['vox_2_sum'] = None
results_df['vox_2_density'] = None
results_df['vox_3_average'] = None
results_df['vox_3_sum'] = None
results_df['vox_3_density'] = None
results_df['full_voxel_dimensions'] = None
results_df['array_side_length_mm'] = None
results_df['vox_resolution_micron'] = None

#Run through logbook and update each row with any supplemental data 1) as available or 2) as called for by user above in 'Main Settings'
for idx, name in enumerate(print_names):
    print('#'*50)
    print(f'Processing logbook row {idx + print_name_start_row}: {name}')

    try:
        #Pull logbook row
        this_logbook_row = logbook_df[logbook_df['Name'] == name]
    
        #Generate a generalized and homogenous structure dict
        diw_structure_dict = structure_dict_from_logbook_row(this_logbook_row)
          # print a quick summary of structure
        print(f"Structure {diw_structure_dict['part_structure']}; Strand {diw_structure_dict['part_skin_nozzle_size']}\{diw_structure_dict['part_layer_nozzle_size']}; Angle {diw_structure_dict['part_angular_offset']}; Offset {diw_structure_dict['part_lateral_offset']}")
        
        ############################################################################################################################
        #####  Add values from opacity/volumetric prediction from logbook structure ################################################
        ############################################################################################################################
        
        #Check if this is a new structure or not; don't update 'unique_structure_dict' as we're just checking
        is_unique_bool,unique_structure_id,_ = match_structure_to_unique_name(diw_structure_dict, unique_structure_dict)
        print(f"\t Structure {diw_structure_dict['part_structure']} flagged as {unique_structure_id}")
        logbook_df.loc[logbook_df['Name']==name, 'unique_structure_id'] = unique_structure_id

        for structure_key in update_key_names:
            #Try and find any of the expected analyses values to add to the logbook row for that structure
            try:
                update_value = waffledata_df[waffledata_df['index']==unique_structure_id][structure_key]
                logbook_df.loc[logbook_df['Name']==name, structure_key] = update_value
            #Make some exceptions for column names that change for clarity in the final 'An
            except KeyError:
                print(f"Key '{unique_structure_id}' does not have sub-key '{structure_key}'")


        ############################################################################################################################
        #####  Add mechanical data updating  #######################################################################################
        ############################################################################################################################

        if update_with_mechanical_data:
            pass

    except Exception as e:
        print()
        print(f"Updating of logbook row with supplemental data failed with error: {e}")
        print()
        print(f"\t {traceback.format_exc()}")

#Save DataFrame to CSV with a similar name
savebook_filename = logbook_filepath.replace("AutomatedAnalysis", "UpdatedResults")
savebook_filepath = os.path.join(os.path.dirname(logbook_filepath), savebook_filename)
with open(savebook_filepath, 'w') as file:
    logbook_df.to_csv(file, index= False, lineterminator='\n')
