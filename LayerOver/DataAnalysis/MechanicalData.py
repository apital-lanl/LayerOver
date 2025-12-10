"""
2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2025-11-25
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Module for pulling and analyzing mechanical data from Excel and CSV files. 

"""

#Import libraries
import csv
import os
import re
from tkinter import filedialog, Tk


  # import other LayerOver modules


#Define variables
example_namerows = {
    'mech-type1': ['Index','Time (sec)', 'Load (kN)', 'Crosshead (mm)', 'PrimaryExtension (mm)', 'Stress (kPa)', 'Strain (mm/mm)', 'Specimen Height (in)'],
    'mech-type2': ['Index','Time (sec)', 'Load (N)', 'Stress (kPa)', 'Gap (mm)', 'GapEx1 (mm)', 'Extension (mm)', 'Extension Ex1 (mm)', 'Strain (mm/mm)', 'Strain Ex1 (mm/mm)']
    }


######################################################################################################################################
##### Generic Utilities  #############################################################################################################
######################################################################################################################################

def walk_directory_for_mech_files(directory=None):
    if not directory:
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()

    elif os.path.isdir(directory):
        #Walk the directory
        for root, dirs, files in os.walk(directory, topdown = True):
            for file in files:
                #Instantialize variables
                this_filepath = os.path.join(root, file)
                immediate_directory_basename = os.path.basename(root)

                if file.endswith('.csv'):
            
                    #Do the filename stuff
                    this_name = file.replace('.csv', '').lower()
                    this_filepath = os.path.join(root, file)
                    csv_filepath_list.append(this_filepath)
                    csv_directory_basenames.append(immediate_directory_basename)
            
                      #try splitting the filename by hyphen or underscore
                    name_split_type = 'unknown'
                    name_split_list_under = this_name.split('_')
                    name_split_list_hyphen = this_name.split('-')
                    if len(name_split_list_under)<2 and len(name_split_list_hyphen)<2:
                        #TODO: add a handler here
                        csv_names_list.append(this_name)
                        name_parser_list.append(name_split_type)
                    elif len(name_split_list_under) < len(name_split_list_hyphen):
                        name_split_type = 'hyphen'
                        this_name = name_split_list_hyphen[0]
                        for part in name_split_list_hyphen[1::]:
                            this_name = str(this_name + str('_'+part))
                        csv_names_list.append(this_name)
                        name_parser_list.append(name_split_type)
                    elif len(name_split_list_under) == len(name_split_list_hyphen):
                        name_split_type = 'ambiguous'
                        csv_names_list.append(this_name)
                        name_parser_list.append(name_split_type)
                    else:
                        name_split_type = 'underscore'
                        csv_names_list.append(this_name)
                        name_parser_list.append(name_split_type)
                
                elif file.endswith('.xlsx'):

    else:
        print("Bad directory; please select an appropriate directory with mechanical data.")
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()


def check_for_namerow(csv_filepath, namerow_example = None):
    result = {
        'likely_name_row': 0,
        'initial_garbage': False,
        'text_rows': {}
    }
    
    with open(csv_filepath, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
        
    max_text_count = 0
    
    for row_index, row in enumerate(rows):
        #Check if each cell in the row contains only numbers
        text_cells = [cell for cell in row if not re.match(r'^-?\d*\.?\d+$', cell.strip())]
        
        #If there are >0 text cells, check if it's a name row or not and what kind
        if text_cells:
            result['text_rows'][row_index] = text_cells

            
            if len(text_cells) > max_text_count:
                max_text_count = len(text_cells)
                result['likely_name_row'] = row_index
    
    result['initial_garbage'] = result['likely_name_row'] > 0
    
    return result


######################################################################################################################################
##### Functions for specific mechanical data  ########################################################################################
######################################################################################################################################

def parse_mech_data_filename(filename):
    '''
    Description: Parse a filename to pull DIW-specific identifiers and names
    INPUT:  
        'filename'      simple filename string; can be whole filepath
    ACTION:
        -find the way the filename is structured
        -split the filename intelligently
        -report results as a dictionary with a standard structure
    OUTPUT:
        'filename_dict' 
    '''
    
    #Initialize variables
    filename_dict = {
        'raw_filename': '',
        'parse_type':'',

        }
    #Make sure the filename is just the basename
    filename = os.path.basename(filename)

    
    return filename_dict