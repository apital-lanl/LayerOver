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
import pandas as pd
import re
from tkinter import filedialog, Tk


  # import other LayerOver modules


#Define variables
example_namerows = {
    'mech-generic': ['Load', 'Stress', 'Strain'],
    'mech-type1': ['Index','Time (sec)', 'Load (kN)', 'Crosshead (mm)', 'PrimaryExtension (mm)', 'Stress (kPa)', 'Strain (mm/mm)', 'Specimen Height (in)'],
    'mech-type2': ['Index','Time (sec)', 'Load (N)', 'Stress (kPa)', 'Gap (mm)', 'GapEx1 (mm)', 'Extension (mm)', 'Extension Ex1 (mm)', 'Strain (mm/mm)', 'Strain Ex1 (mm/mm)']
    }

blank_mech_file_metadata_dict = {
    'filename':'',
    'filepath':'',
    'filetype':'',
    'immediate_parent_directory':'',
    'parse_type':'',
    'clean_filename':'',
    'root_print_name?':False,
    'print_name':'',
    'print_project_name':'',
    'thickness': 0,
    'density': 0,
    'material': '',

    }


######################################################################################################################################
##### Generic Utilities  #############################################################################################################
######################################################################################################################################

def process_directory_for_mech_files(directory=None):
    #Initialize variables
    mech_data_dict = {}
      # if no directory is added as input, select one
    if not directory:
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()


def walk_directory_for_mech_files(directory=None):
    '''
    Description: Capture all CSV and XLSX files within a directory. Parse the spreadsheet and return summary data and metadata.
    '''

    #Initialize variables
    mech_filedata_dict = {}
      # if no directory is added as input, select one
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
                this_dict = blank_mech_file_metadata_dict.copy()

                if file.endswith('.csv'):
                    #Do the filename stuff
                    this_dict['filetype'] = '.csv'
                    this_dict['filename'] = file
                    this_filepath = os.path.join(root, file)
                    this_dict['filepath'] = this_filepath
                    this_dict['immediate_parent_directory'] = immediate_directory_basename
            
                      #try splitting the filename by hyphen or underscore
                    parse_dict = parse_mech_data_filename(file)  #returns dict with 'raw_filename', 'parse_type', 'clean_filename'
                    this_dict['parse_type'] = parse_dict['parse_type']
                    this_dict['clean_filename'] = parse_dict['clean_filename']

                    #Check CSV contents

                
                elif file.endswith('.xlsx'):
                    #Do the filename stuff
                    this_dict['filetype'] = '.xlsx'
                    this_dict['filename'] = file
                    this_filepath = os.path.join(root, file)
                    this_dict['filepath'] = this_filepath
                    this_dict['immediate_parent_directory'] = immediate_directory_basename
            
                      #try splitting the filename by hyphen or underscore
                    parse_dict = parse_mech_data_filename(file)  #returns dict with 'raw_filename', 'parse_type', 'clean_filename'
                    this_dict['parse_type'] = parse_dict['parse_type']
                    this_dict['clean_filename'] = parse_dict['clean_filename']
    else:
        print("Bad directory; please select an appropriate directory with mechanical data.")
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()

    return mech_filedata_dict


def check_for_namerow(spreadsheet_filepath, 
                      namerow_example = example_namerows,
                      row_limit = 50):
    '''
    Description: Generic function to open a spreadsheet (CSV or XLSX), find likely namerow based on an example or largest text-containing row, 
        and return a dict with metadata from the namerow search.

    INPUT:
        'spreadsheet_filepath'  str; filepath for a CSV or XLSX file. If XLSX, only one sheet will be returned.
      (optional)
        'namerow_example'       dict; each key is a type of file, and the list stored at that key is an example namerow
        'row_limit'             int; maximum number of rows to look for column header text in
    ACTIONS:
        - 
    '''
    
    #Initialize variables
    namerow_dict = {
        'likely_name_row': 0,
        'initial_garbage': False,
        'data_type_guess': '',
        'data_type_match_count': 0,
        'text_rows': {}
        }
      # create a dict to hold matches for each example type
    max_example_matches = {}
    for key in list(namerow_example.keys()):
        if 'generic' in key:
            generic_key = key
        max_example_matches[key] = 0
    min_text_count = 0   #minimum number of text rows required to be a namespace
    ideal_text_count = 0   #likely number of namerows based on number of columns is data rows
    
    if spreadsheet_filepath.endswith('.csv'):
        with open(spreadsheet_filepath, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)

    elif spreadsheet_filepath.endswith('.xlsx'):
        pass
        
    #Not currently used; Hard-coded settings for match quality
    max_row_index_to_consider = row_limit   #assume any rows below this can't possibly have column labels
    max_match_similarity = 0
        
    #Consider each row for text
    for row_index, row in enumerate(rows):
        #Check if each cell in the row contains only numbers
        text_cells = [cell for cell in row if not re.match(r'^-?\d*\.?\d+$', cell.strip())]
        
        #If there are >0 text cells, check if it's a name row or not, and what kind
        if text_cells and (row_index < max_row_index_to_consider):
            #Initialize new variables
            namerow_dict['text_rows'][row_index] = {
                'text_cels': text_cells,
                'example_similarities': {} 
                }
            match_type_guess = 'unknown'
            max_similarity_match = 0
                # make everything lowercase
            text_cells = [text.lower() for text in text_cells]

            #Pull each example and get a similarity match
            for this_key in list(namerow_example.keys()):
                this_match_count = 0
                these_examples = namerow_example[this_key]
                #Nested For loops to allow for fuzzy matching and because these are tiny datasets that are only called sparingly
                #TODO: add fuzzy matching
                for item in these_examples:
                    item = item.lower()
                    for text in text_cells:
                        if item in text:
                            this_match_count += 1
                #Add the results to the growing dict
                namerow_dict['text_rows'][row_index]['example_similarities'][this_key] = this_match_count
                if this_match_count > max_similarity_match:
                    max_similarity_match = this_match_count

            
        #TODO: add possible flag and recognition for non-text rows (to aid in datatyping later)
        else:
            pass

    for text_row_idx in list(namerow_dict['text_rows'].keys()):
        this_row_dict = namerow_dict['text_rows'][text_row_idx]['example_similarities']
        for example_type in list(max_example_matches.keys()):
            if this_row_dict[example_type] > max_example_matches[example_type]:
                max_example_matches[example_type] = this_row_dict[example_type]
                if this_row_dict[example_type] > max_match_similarity:
                    namerow_dict['data_type_guess'] = example_type
                    namerow_dict['likely_name_row'] = text_row_idx
                    max_match_similarity = this_row_dict[example_type]
           
        
    

    
    #Set a boolean flag for initial garbage
    #TODO: add characterization of the garbage
    namerow_dict['initial_garbage'] = namerow_dict['likely_name_row'] > 0

    return namerow_dict


def pull_last_mechanical_replicate(data_df, data_dict = None):
    '''
    Description: Take cyclic load test data and just return the final loading/unloading cycle as separate columns.
    INPUT:
        'data_df'       pandas DataFrame; Contains one 'Stress' and one 'Strain' column
                        if not a dataframe, raise error
                        if more than one stress or strain column, raise error
    ACTION:
        -lorem
    OUTPUT:
        'mech_df'       pandas DataFrame; 'Index', 'Stress', 'Strain (loading)', 'Strain (unloading)

    '''
    #Initialize variables
    mech_dict = {
        'testing_index': [],
        'stress_data': [],
        'strain_data_loading': [],
        'strain_data_unloading': []
        }
    mech_units = {
        'testing_index': 'int'
        'stress_data': 'kN',
        'strain_data_loading': r'mm/mm'
        'strain_data_unloading': r'mm/mm'
        }

    return mech_df


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
        'clean_filename':''
        }
      # make sure the filename is just the basename
    filename = os.path.basename(filename)
    filename_dict['raw_filename'] = filename
      # split out filename from filetype
    filename = os.path.splitext(filename)[0]
      # initialize differently parsed names
    name_split_list_under = this_name.split('_')
    name_split_list_hyphen = this_name.split('-')
    if len(name_split_list_under)==0 and len(name_split_list_hyphen)==0:
        #If for some reason the input string is faulty, report as a failure; should never happen
        filename_dict['parse_type'] = 'failure'
        filename_dict['clean_filename'] = filename

    elif len(name_split_list_under)<2 and len(name_split_list_hyphen)<2:
        #TODO: add a handler here
        #If both underscore and hyphen parsing have only 1 array element, the filename is a single word or character string
        filename_dict['parse_type'] = 'direct'
        filename_dict['clean_filename'] = filename

    elif len(name_split_list_under) < len(name_split_list_hyphen):
        #At least two parts of the filename are hyphenated and fewer are underscored
        filename_dict['parse_type'] = 'hyphen'
          # turn hyphenated into underscored to match logbook
        this_name = name_split_list_hyphen[0]
        for part in name_split_list_hyphen[1::]:
            this_name = str(this_name + str('_'+part))
        filename_dict['clean_filename'] = this_name
    
    elif len(name_split_list_under) > len(name_split_list_hyphen):
        #Should be 'standard' underscored filename
        filename_dict['parse_type'] = 'underscore'
        filename_dict['clean_filename'] = filename    #matches logbook format, just pass filename through

    else:
        #Should be very rare case that there are >2 underscores AND hyphens or other parsing error not captured above
        filename_dict['parse_type'] = 'ambiguous'
          # turn hyphenated into underscored to match logbook
        filename_dict['clean_filename'] = filename

    return filename_dict