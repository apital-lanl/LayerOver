"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   20YY-MM-DD
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Small script utilizing LayerOver modules to parse and open mechanical data files.

"""

#Import outside libraries
from tkinter import Tk, filedialog
import os
import pandas as pd
import matplotlib.pyplot as plt
import traceback
  
#Imports from LayerOver 
from LayerOver.Analysis import MechanicalData as mech


#Specify some display options
show_each_file_results = True               #show outputs from each file that's read
show_displacement_peak_locations = False     #show where mechanical replicate indices have been chosen
displacement_type = 'extension'             #'extension'    'strain'

#Select a file(s) to import
root = Tk()
filenames = filedialog.askopenfilenames(title= "Select mech data files to import", filetypes = [('Spredsheet files', '*.csv *.xlsx')])
root.destroy()
#Select a logbook to import
root = Tk()
directory = filedialog.askdirectory(title= "Select a folder with an updated logbook for comparison")
root.destroy()

#Try and find the 'latest' logbook
    # brute-force a logbook search for a couple of levels of parent directory if none is found in analysis directory
logbook_filepath = mech.get_latest_logbook(directory)   #Returns a str if successful and NONE if not successful
    # immediate parent directory
if not logbook_filepath:
    parent_directory = os.path.dirname(os.path.dirname(directory))
    logbook_filepath = mech.get_latest_logbook(parent_directory)
    # grand-parent directory
if not logbook_filepath:
    parent_directory = os.path.dirname(os.path.dirname(parent_directory))
    logbook_filepath = mech.get_latest_logbook(parent_directory)
    # open the logbook as a pandas DataFrame

#Try and open the latest filepath
if logbook_filepath:
    logbook_df = mech.open_logbook(logbook_filepath,
                                target_excel_sheetname = None)
else:
    #TODO: need to add final recourse here if no logbook has been found
    print("Can't find a logbook to open.")

#Run each selected file
for filename in filenames:
    print()
    print("#"*50)
    print(os.path.basename(filename))
    
    try:
        this_file_dict = mech.process_mechanical_file_against_logbook(filename, logbook_df,
                                                                alternate_save_directory = directory,
                                                                show_displacement_peaks = show_displacement_peak_locations,
                                                                show_each_file_results = show_each_file_results,
                                                                save_last_replicate_graph = show_each_file_results,
                                                                displacement_column = displacement_type)

        # 'this_file_dict' keys and values:
        #-----------------------------------------------------------
        # 'filename'                    str
        # 'filepath'                    str
        # 'filetype'                    str; 'csv' or 'xlsx'
        # 'immediate_parent_directory'  str
        # 'parse_type'                  str; 'hyphen', 'underscore', 'UNK'
        # 'clean_filename'              str
        # 'data_namerow_dict'           dict; info on which (if any) rows contain likely column names for mechanical data (i.e. 'stress', 'strain')
        # 'pandas_readable'             bool; did parsing into a pd.DataFrame work?
        # write results to the output row


    except Exception as le:
        #'le' stands for 'loop error'
        print(f"Exception thrown: {le}")
        print(traceback.format_exc())
        print()
        print("Failure to parse; moving on to next file.")
        print("_"*50)
        print()