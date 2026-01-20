"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
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

TODO:
    -Need to implement XLSX parsing. Sticking point is going through each spreadsheet and I don't want to deal with it right now.

"""

#Import libraries
import csv
import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from tkinter import filedialog, Tk
from scipy.signal import find_peaks

  # import other LayerOver modules or module parts
from LayerOver.PSPP.DIWStructure import blank_diw_structure_dict

#Define variables
#Define hard-coded thresholds and setting values
default_stress_threshold = 0.2   #in kPa; for silicone elastomers, but should be relatively general
strain_minimum_mask_threshold = -0.1  #minimum strain to accept (<0 to allow for noise at 0 strain)

#Example column names for various types of report from mechanical testing instruments
# used to guess which 1) type of mech data is being parsed, 2) which row in the spreadsheet contains the column names, 
# and 3) which columns to pull for 'stress' and 'strain' data.
example_namerows = {
    'mech-generic': ['Load', 'Stress', 'Strain'],
    'mech-type1': ['Index','Time (sec)', 'Load (kN)', 'Crosshead (mm)', 'PrimaryExtension (mm)', 'Stress (kPa)', 'Strain (mm/mm)', 'Specimen Height (in)'],
    'mech-type2': ['Index','Time (sec)', 'Load (N)', 'Stress (kPa)', 'Gap (mm)', 'GapEx1 (mm)', 'Extension (mm)', 'Extension Ex1 (mm)', 'Strain (mm/mm)', 'Strain Ex1 (mm/mm)']
    }

#Details for parsing files of the types in 'example_namerows'
#'__exclusion_terms' are lists of terms that disqualify a name. I.E. 'Strain Ex1 (mm/mm)' is a valid 'strain' name for 'mech-type2', but we only want the 'Strain (mm/mm)' column and 
#       and should ignore any 'strain' column with ' Ex1' also in the name.
namerow_example_parsing_dict = {
    'mech-generic': {
        'stress_columnn_name': 'Stress',
        'strain_column_name': 'Strain',
        'stress_exclusion_terms': ['peak'],
        'strain_exclusion_terms': [' Ex1 ']
        },
    'mech-type1': {
        'stress_columnn_name': 'Stress (kPa)',
        'strain_column_name': r'Strain (mm/mm)',
        'stress_exclusion_terms': ['peak'],
        'strain_exclusion_terms': []
        },
    'mech-type2': {
        'stress_columnn_name': 'Stress (kPa)',
        'strain_column_name': r'Strain (mm/mm)',
        'stress_exclusion_terms': ['peak'],
        'strain_exclusion_terms': []
        }
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

default_units_dict = {
    'mm': {
        'return':'mm',
        'exact_matcth': '(mm)',
        'alternates':[r'mm/mm']
        },
    r'mm/mm': {
        'return':r'mm/mm',
        'exact_matcth': r'(mm/mm)',
        'alternates':[]
        },
    'N': {
        'return':'N',
        'exact_match':'(N)',
        'alternates':['kN']
        },
    'kN': {
        'return':'kN',
        'exact_match':'(kN)',
        'alternates':[]
        },
    'Pa': {
        'return':'Pa',
        'exact_match':'(Pa)',
        'alternates':['kPa']
        },
    'kPa': {
        'return':'kPa',
        'exact_match':'(kPa)',
        'alternates':[]
        },
    'sec': {
        'return':'sec',
        'exact_match':'(sec)',
        'alternates':[]
        },
    'in': {
        'return':'in',
        'exact_match': '(in)',
        'alternates':[]
        },
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

    #Get all CSVs in the directory
    mech_filedata_dict = walk_directory_for_mech_files(directory)
    mech_file_keys = list(mech_filedata_dict.keys())

    #Open each file and process
    for filename in mech_file_keys:
        
        #Pull the values for this particular mechanical data file
        this_filedata_dict = mech_filedata_dict[filename]
        #'this_filedata_dict' values: 
        #   'filename'                      str         populated
        #   'filepath'                      str         populated
        #   'filetype'                      str         populated
        #   'immediate_parent_directory'    str         populated
        #   'parse_type'                    str         populated
        #   'clean_filename'                str         populated
        #   'root_print_name?':False,       bool        defualt
        #   'print_name':'',                str         defualt
        #   'print_project_name':'',        str         defualt
        #   'thickness': 0,                 int         default
        #   'density': 0,                   int/float   default
        #   'material': '',                 str         defualt
        this_filepath = this_filedata_dict['filepath']
        this_filetype = this_filedata_dict['filetype']
        this_mechdata_parsetype = this_filedata_dict['filetype']   #what kind of instrument ouput (i.e. column names) should the file be opened with?
        
        #Try to open the file, get the data, and compare to 
        data_dict = parse_mech_data_fromcsv(this_filepath)
        data_df = data_dict['dataframe']

        logbook_dict = check_logbook_for_printname(printname)


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

                    ################################################
                    # TODO: add some parsing flags to 'this_dict'
                    ################################################
                    
                    #Add dict to growing return dict
                    mech_filedata_dict.update({file: this_dict})


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

                    ################################################
                    # TODO: add some parsing flags to 'this_dict'
                    ################################################

                    #Add dict to growing return dict
                    mech_filedata_dict.update({file: this_dict})
    else:
        print("Bad directory; please select an appropriate directory with mechanical data.")
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()

    return mech_filedata_dict


def check_for_namerow(spreadsheet_filepath, 
                      namerow_example = example_namerows,
                      row_limit = 50,
                      show_peaks = True):
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
        'successful_parse': False,
        'likely_name_row': 0,
        'initial_garbage': False,
        'data_type_guess': '',
        'data_type_match_count': 0,
        'filetype': '',
        'data_sheetname':'',
        'text_rows': {}
        }
      # create a dict to hold matches for each example type
    max_example_matches = {}
    for key in list(namerow_example.keys()):
        if 'generic' in key:
            #TODO: add a comparison between 'generic' matches and specific match types
            generic_key = key
        max_example_matches[key] = 0
    min_text_count = 0   #minimum number of text rows required to be a namespace
    ideal_text_count = 0   #likely number of namerows based on number of columns is data rows
    
    if spreadsheet_filepath.endswith('.csv'):
        with open(spreadsheet_filepath, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)
        namerow_dict['filetype'] = 'csv'

    elif spreadsheet_filepath.endswith('.xlsx'):
        #Open excel file as a pandas object
        namerow_dict['filetype'] = 'xlsx'
        try:
            excel_file = pd.ExcelFile(spreadsheet_filepath)
        except:
            print("Failure to read Excel file: this usually means it's corrupted.")
            return namerow_dict

        #If only one sheet, return the rows from that sheet.
        if len(excel_file.sheet_names) ==1:
            sheet_name = list(excel_file.sheet_names)[0]
            excel_df = pd.read_excel(excel_file, sheet_name = sheet_name, nrows=10)
            rows = excel_df.values.tolist()
        #If more than one sheet, try and find the best one
        elif len(excel_file.sheet_names) > 1:
            for idx, sheet_name in enumerate(excel_file.sheet_names):
                this_df = pd.read_excel(excel_file, sheet_name = sheet_name)
                these_rows = this_df.values.tolist()
                cur_max_match_score = 0
                default_namerow = example_namerows['mech-generic']
                namerow_max_score = len(default_namerow)
                for row_idx, row in enumerate(these_rows):
                    text_cells = [cell for cell in row if not re.match(r'^-?\d*\.?\d+$', cell.strip())]
                    this_row_score = 0
                    if len(text_cells) > 0:
                        for this_cell in text_cells:
                            this_check = sum([1 for name in default_namerow if this_cell.lower() in name.lower()])
                            this_row_score += this_check
                    if this_row_score > cur_max_match_score:
                        cur_max_match_score = this_row_score
                #If every default name gets at least 1 match, assume this is the sheet
                #NOTE: if more than 1 sheet has mech data, this will return the last mech data sheet's info
                if cur_max_match_score >= namerow_max_score:
                    namerow_dict['data_sheetname'] = sheet_name
                    rows = these_rows  #pass to loop below
                elif cur_max_match_score >0:
                    #TODO: if some matches are encountered, handle to find best sheet
                    pass

        else:
            print()
            print(f"Failure to parse {spreadsheet_filepath}")
            return namerow_dict
        
    #Not currently used; Hard-coded settings for match quality
    max_row_index_to_consider = row_limit   #assume any rows below this can't possibly have column labels
    max_match_similarity = 0
        
    #Consider each row for text
    for row_index, row in enumerate(rows):
        #Check if each cell in the row contains only numbers
        text_cells = [cell for cell in row if not re.match(r'^-?\d*\.?\d+$', str(cell).strip())]
        
        #If there are >0 text cells, check if it's a name row or not, and what kind
        if text_cells and (row_index < max_row_index_to_consider):
            #Initialize new variables
            namerow_dict['successful_parse'] = True
            namerow_dict['text_rows'][row_index] = {
                'text_cels': text_cells,
                'example_similarities': {} 
                }
            match_type_guess = 'unknown'
            max_similarity_match = 0
                # make everything lowercase
            text_cells = [str(text).lower() for text in text_cells]

            #Pull each example and get a similarity match
            for this_key in list(namerow_example.keys()):
                this_match_count = 0
                these_examples = namerow_example[this_key]
                #Nested For loops to allow for fuzzy matching and because these are tiny datasets that are only called sparingly
                #TODO: add fuzzy matching; currently only matches exactly (including spacing) but ignores case (i.e. upper/lower)
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

    #Find the best guess for mechanical data format and 
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


def pull_mechanical_replicates(data_df, data_dict = None,
                                   show_peaks = True,
                                   stress_threshold = None,
                                   strain_zero_offset = 10,
                                   strain_min_thresh = None,
                                   report_nonnegative_strain = False):
    '''
    Description: Take cyclic load test data and just return the final loading/unloading cycle as separate columns.
    INPUT:
        'data_df'       
            pandas DataFrame; Contains one 'Stress' and one 'Strain' column
            if not a dataframe, raise error
            if more than one stress or strain column, raise error

        (OPTIONAL)
        'stress_threshold'
            Stress value for a moving average at which strain is evaluated as starting to affect material.
            Strain 0 is reset to a value a little before this threshold is met based on 'strain_zero_offset'.
        'strain_zero_offset'
            Number of datapoints before stress is starting to register that is reset as the new '0 strain'
            Set by first replicate and all future replicates take this strain value.
                        
    ACTION:
        -lorem
    OUTPUT:
        'mech_df'       pandas DataFrame; 'Index', 'Stress', 'Strain (loading)', 'Strain (unloading)

    Notes:
        -Typical strain data is ~1-5 microns (.001-.005 mm) difference between prior point and next point; any less than ~1e-4 change is probrably an inflection point

    '''
    #Initialize variables
    replicate_dict = {
        'replicate_parse_success': False,
        'number_of_replicates': 0,
        'last_strain_peak_index': 0,
        'last_strain_valley_index': 0,
        'peak_to_valley_index_diff': 0,
        'replicate_data': {
            }
        }
    mech_dict = {
        'testing_index': [],
        'all_strain_data': [],
        'all_stress_data': [],
        }

    data_dict = {
        'mech_units_info': {
            'testing_index': 'UNK',
            'stress_units': 'UNK',
            'strain_units': 'UNK'
            },
        }
    columns = list(data_df.columns)
    stress_flag = False
    strain_flag = False

    #Check for 'data_df' data
    #TODO: add DataFrame type check

    #Flip and coerce strain values 

    #Try pulling stress and strain data
    for column_name in columns:
        try:
            unit_guess = try_and_guess_units(column_name)
        except:
            print()
            print(f"Failure to guess units on column name: {column_name}")
            unit_guess = 'UNK'
        
        if ('stress' in column_name.lower()) and not stress_flag:
            mech_dict['all_stress_data'] = data_df[column_name]
            print()
            print(f"Added {column_name} as stress data")
            stress_flag = True
            data_dict['mech_units_info']['stress_units'] = unit_guess
            stress_units = unit_guess
            stress_col_name = f"Stress ({stress_units})"
        elif ('stress' in column_name.lower()):
            print
            print(f"Multiple stress columns passed. Ignoring {column_name}")
        
        if ('strain' in column_name.lower()) and not strain_flag:
            mech_dict['all_strain_data'] = data_df[column_name]
            print()
            print(f"Added {column_name} as strain data")
            strain_flag = True
            data_dict['mech_units_info']['strain_units'] = unit_guess
            strain_units = unit_guess
            strain_col_name = f"Strain ({strain_units})"
        elif ('strain' in column_name.lower()):
            print
            print(f"Multiple strain columns passed. Ignoring {column_name}")

    #Use peak finding to grab last loading/unloading cycl
    raw_df = pd.DataFrame(data= {stress_col_name: mech_dict['all_stress_data'],
                                  strain_col_name: mech_dict['all_strain_data']})
    
      # find all peaks in strain data
    peaks, _ = find_peaks(raw_df[f'Strain ({strain_units})'].values)
    replicate_dict['number_of_replicates'] = len(peaks)
    
    if len(peaks) == 0:
        print("No peaks found in strain data")
    #Use last peak to define variables
    elif len(peaks) <2:
        #Find the last peak index
        last_peak_index = peaks[-1]
          # assume last peak index is followed by a 'standard' strain region, so the straining cycle indexes can be backed-out 
        cycle_index_diff_guess = raw_df[f'Strain ({strain_units})'].shape[0] - last_peak_index
        last_valley_guess_index = last_peak_index - cycle_index_diff_guess
        #Assign values to the return dict
        replicate_dict['replicate_parse_success'] = True
        replicate_dict['last_strain_peak_index'] = last_peak_index
        replicate_dict['last_strain_valley_index'] = last_valley_guess_index
        replicate_dict['peak_to_valley_index_diff'] = cycle_index_diff_guess
    else:
        #Find the last peak index
        last_peak_index = peaks[-1]
        cycle_index_diff_guess = (peaks[-1]-peaks[-2])//2  #get difference between last two peaks
        last_valley_guess_index = last_peak_index - cycle_index_diff_guess
        #Assign values to the return dict
        replicate_dict['replicate_parse_success'] = True
        replicate_dict['last_strain_peak_index'] = last_peak_index
        replicate_dict['last_strain_valley_index'] = last_valley_guess_index
        replicate_dict['peak_to_valley_index_diff'] = cycle_index_diff_guess
    
    #Show the peak locations if prompted and replicate parsing was successful
    if show_peaks and replicate_dict['replicate_parse_success']:
        #Pull relevant indices and graph values for replicate indices
        stress_max = mech_dict['all_stress_data'].max()
        strain_max = mech_dict['all_strain_data'].max()

        #Generate the graph
        plt.figure(figsize = (10,10))
        plt.scatter(list(range(mech_dict['all_strain_data'].shape[0])), mech_dict['all_strain_data'])
        plt.title(f"Strain peaks found")
        for peak in peaks:
            plt.scatter(peak, mech_dict['all_strain_data'][peak], marker = 'x', s = 200, color = 'gray')
          # start of last loading replicate
        plt.plot([last_valley_guess_index, last_valley_guess_index], [0, strain_max], color='r', linewidth=3, alpha=0.6)
          # start of last unloading replicate
        plt.plot([last_peak_index, last_peak_index], [0, strain_max], color='g', linewidth=3, alpha=0.6)
        plt.show()
    #If replicate parsing not successful, show a graph that might hint at why
    elif show_peaks:
        plt.figure(figsize = (10,10))
        plt.scatter(list(range(mech_dict['all_strain_data'].shape[0])), mech_dict['all_strain_data'])
        plt.title(f"Failed to find strain peaks")
        plt.show()
    
    #Loading cycles start at (peak_idx-valley_idx), Unloading cycles end at (peak_idx+valley_idx)
    for replicate_idx, peak_idx in enumerate(peaks):
        load_start_idx = peak_idx-cycle_index_diff_guess
        unload_end_idx = peak_idx+cycle_index_diff_guess
        #Make sure indices are in-bounds for data size
        if load_start_idx<0:
            load_start_idx = 0
          # should be the same size, but just in case take the minimum of the two columns dimensions
        smallest_data_index = min(raw_df[strain_col_name].shape[0], raw_df[stress_col_name].shape[0])-1  #subtract 1 for 0-start indexing
        if unload_end_idx > smallest_data_index:
            unload_end_idx = smallest_data_index

        #Find strain minimum where stress actually starts increasing above a threshold
        if replicate_idx == 0:
              # if no setting is passed, use module default
            if not stress_threshold:
                stress_threshold = default_stress_threshold

            #Grap the loading curve for analysis of strain 0 offset
            strain_loading = raw_df[strain_col_name].iloc[load_start_idx:peak_idx]
            stress_loading = raw_df[stress_col_name].iloc[load_start_idx:peak_idx]

            #We'll only consider the first loading cycle
            stress_loading_avg = stress_loading.rolling(3, center=True, min_periods = 1).mean()
            valid_stress_mask = stress_loading_avg >= stress_threshold
            first_valid_stress_index = valid_stress_mask.idxmax()
            first_valid_strain_index = first_valid_stress_index-strain_zero_offset
            if first_valid_strain_index < 0:
                first_valid_strain_index = 0
            #Get the strain value at new predicted '0 strain' value
            first_valid_strain = raw_df[stress_col_name].iloc[first_valid_strain_index]
            if first_valid_strain <0:
                first_valid_strain = 0
            #Reset the entire strain data
            raw_df[strain_col_name] = raw_df[strain_col_name]- first_valid_strain

        #Grab cycle data
        strain_loading = raw_df[strain_col_name].iloc[load_start_idx:peak_idx]
        stress_loading = raw_df[stress_col_name].iloc[load_start_idx:peak_idx]
        strain_unloading = raw_df[strain_col_name].iloc[peak_idx:unload_end_idx]
        stress_unloading = raw_df[stress_col_name].iloc[peak_idx:unload_end_idx]
        if report_nonnegative_strain:
            #If no min strain threshold is passed, use module default
            if not strain_min_thresh:
                strain_min_thresh = strain_minimum_mask_threshold
            #Get pandas Series mask for values of strain above threshold
            loading_mask = strain_loading[strain_loading>= strain_min_thresh]
            unloading_mask = strain_unloading[strain_unloading>= strain_min_thresh]
            #Apply non-negative mask to get appropriate values only
            strain_loading = strain_loading[loading_mask]
            stress_loading = stress_loading[loading_mask]
            strain_unloading = strain_unloading[unloading_mask]
            stress_unloading = stress_unloading[unloading_mask]
    
        #Create output DataFrame for this replicate
        output_df = pd.DataFrame({
            'strain_data_loading': strain_loading,
            'strain_data_unloading': strain_unloading,
            'stress_data_loading': stress_loading,
            'stress_data_unloading': stress_unloading
            })
    
        #TODO: fix this; not sure why it kept throwing errors and it doesn't really matter so I've moved on  
          # reindex to ensure all columns have the same length
        max_length = max(output_df[col].shape[0] for col in output_df.columns)
        #output_df = output_df.reindex(range(max_length))

        #Store replicate data with replicate index (starts at 1)
        replicate_dict['replicate_data'][replicate_idx+1] = output_df

    return replicate_dict


######################################################################################################################################
##### Functions for specific mechanical data  ########################################################################################
######################################################################################################################################

def parse_mech_data_filename(filename):
    '''
    Description: Parse a filename to pull DIW-specific identifiers and names.
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

    #Try different name parsing schemes
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


def parse_mech_data_fromcsv(mech_data_filepath,
                            data_dict = None,
                            special_example_namerow_dict = None):
    '''
    Description:
        Open a CSV file, check for stress/strain data, and output a pandas DataFrame if that data exists.
        Also homogenizes the strain data prior to exporting.

    INPUT:
        (optional)
        special_example_namerow_dict    dict; at least one key:list pair with list of column name strings in a special mechanical data file example
        data_dict                       dict; not used for CSV parsing currently, but added as future option and for consistency with XLSX parsing

    '''

    #Initialize variables
    file_output_dict = {
        'parse_status': {
            'file_pandas_readable': False,
            'namerow_good': False,
            'namerow_idx': 0,

            },
        'filename': os.path.basename(mech_data_filepath),
        'filepath': mech_data_filepath,
        'dataframe': None,
        }

    #Try to open the file, find the namerow, and match to expected mechanical data formats
      # if a namerow example has been passed, use that. Otherwise use the defualt at the header of this file.
    if not special_example_namerow_dict:
        namerow_dict = example_namerows
        
    else:
        namerow_dict = special_example_namerow_dict
      # check for a row with mechanical data column names from the example dict
    namerow_dict = check_for_namerow(mech_data_filepath, 
                                          namerow_example = namerow_dict)
    #'namerow_dict' keys:    successful_parse (bool);   likely_name_row (int);   initial_garbage (bool)
    #                        data_type_guess (str);  data_type_match_count (int);  text_rows (dict)
    
    #Try to open the file if namerow check hasn't completely failed
    if namerow_dict['successful_parse'] == True:
        start_row_idx = namerow_dict['likely_name_row']
        mech_data_archetype_guess = namerow_dict['data_type_guess']
        #Try to get parsing details for column names from 'namerow_example_parsing_dict' above; default to 'mech-generic' if things get weird
        try:
            this_parsing_dict = namerow_example_parsing_dict[mech_data_archetype_guess]
        except KeyError:
            this_parsing_dict = namerow_example_parsing_dict['mech-generic']
        stress_name = this_parsing_dict['stress_columnn_name']
        strain_name = this_parsing_dict['strain_column_name']
        stress_name_exclusions = this_parsing_dict['stress_exclusion_terms']
        strain_name_exclusions = this_parsing_dict['strain_exclusion_terms']

        #Handle CSV
        if mech_data_filepath.lower().endswith('.csv'):
            #Open the CSV and find the column name row
            with open(mech_data_filepath, 'r') as file:
                reader = csv.reader(file)
                for row_num, row in enumerate(reader):
                    if row_num == start_row_idx:
                        header_row = row
            
            #Walk through header row and find which column index has the stress and strain columns
            for col_idx, cell_text in enumerate(header_row):
                if stress_name in cell_text:
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in stress_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        stress_col_idx = col_idx
                        stress_col_name = cell_text
                if strain_name in cell_text:
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in strain_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        strain_col_idx = col_idx
                        strain_col_name = cell_text

            #Finally, open the data
            data_df = pd.read_csv(mech_data_filepath, usecols=[stress_col_idx, strain_col_idx], skiprows=(start_row_idx) )

            #Make sure the strain is right-side up (i.e. positive values only), ignore lead-in if not at 0, and set minimum at 0
            try:
                strain_series = data_df[strain_col_name].copy()
            except KeyError:
                #if keyerror, parsing has failed and take a look at the data_df to see what's going on
                print(data_df.head(10))
            
              # get info about strain
            pos_strain_sum = len(strain_series[strain_series>0])
            neg_strain_sum = len(strain_series[strain_series<0])
            this_min_strain = strain_series.min()
            
              # flip strain if negative
            if neg_strain_sum > pos_strain_sum:
                # adj_strain_series = strain_series -this_min_strain
                # adj_strain_series = abs(this_min_strain)-adj_strain_series
                adj_strain_series = strain_series * -1
                strain_series = adj_strain_series  #reset series
            
              # check if a lead-in is artificially skewing 0
            data_start_idx = round(len(strain_series)*0.25)  #ignore first part
            new_min = strain_series[data_start_idx::].min()  #find real '0' without lead-in garbage
            strain_series = strain_series + abs(new_min)  #reset series minimum to a more useful 0

              # re-assign strain series
            data_df[strain_col_name] = strain_series

            #Adjust for any 0-stress lead-in
  

            #Assign everything to the output dict
            file_output_dict['dataframe'] = data_df
            file_output_dict['parse_status']['file_pandas_readable'] = True
            file_output_dict['parse_status']['namerow_good'] = True
            file_output_dict['parse_status']['namerow_idx'] = start_row_idx

        #Provide cross-compatibility just in case
        if mech_data_filepath.lower().endswith('.xlsx'):
            file_output_dict = parse_mech_data_fromxlsx(mech_data_filepath,
                                                        data_dict = data_dict)

    return file_output_dict


def parse_mech_data_fromxlsx(mech_data_filepath,
                            data_dict = None,
                            special_example_namerow_dict = None):
    '''
    Description:
        Open an XLSX file, check for stress/strain data, and output a pandas DataFrame if that data exists.
        Also homogenizes the strain data prior to exporting.

    INPUT:
        (optional)
        special_example_namerow_dict    dict; at least one key:list pair with list of column name strings in a special mechanical data file example
        data_dict                       dict; contains data on which sheet to pull XLSX data from 

    '''

    #Initialize variables
    file_output_dict = {
        'parse_status': {
            'file_pandas_readable': False,
            'namerow_good': False,
            'namerow_idx': 0,

            },
        'filename': os.path.basename(mech_data_filepath),
        'filepath': mech_data_filepath,
        'dataframe': None,
        }

    #Try to open the file, find the namerow, and match to expected mechanical data formats
      # if a namerow example has been passed, use that. Otherwise use the defualt at the header of this file.
    if not special_example_namerow_dict:
        namerow_dict = example_namerows
    else:
        namerow_dict = special_example_namerow_dict
    
    #Check if passed data_dict exists and contains good info on proper sheet to pull. Otherwise just reprocess the file.
    if not data_dict:
        #Dummy check for finding a 'good_sheet_guess' in case previous
        namerow_dict = check_for_namerow(mech_data_filepath, 
                                          namerow_example = example_namerows,
                                          row_limit = 50,
                                          show_peaks = True)
        good_sheetname = namerow_dict['data_sheetname']
    else:
        try:
            namerow_dict = data_dict.copy()
            good_sheetname = data_dict['data_sheetname']
        except:
            print(r'\n', "Cannot parse 'data_dict' input; reprocessing file for good sheet location.")
            namerow_dict = check_for_namerow(mech_data_filepath, 
                                    namerow_example = example_namerows,
                                    row_limit = 50,
                                    show_peaks = True)
            good_sheetname = namerow_dict['data_sheetname']
        

      # check for a row with mechanical data column names from the example dict
    namerow_dict = check_for_namerow(mech_data_filepath, 
                                    namerow_example = namerow_dict)
    #'namerow_dict' keys:    successful_parse (bool);   likely_name_row (int);   initial_garbage (bool)
    #                        data_type_guess (str);  data_type_match_count (int);  text_rows (dict)
    
    #Try to open the file if namerow check hasn't completely failed
    if namerow_dict['successful_parse'] == True:
        start_row_idx = namerow_dict['likely_name_row']
        mech_data_archetype_guess = namerow_dict['data_type_guess']
        #Try to get parsing details for column names from 'namerow_example_parsing_dict' above; default to 'mech-generic' if things get weird
        try:
            this_parsing_dict = namerow_example_parsing_dict[mech_data_archetype_guess]
        except KeyError:
            this_parsing_dict = namerow_example_parsing_dict['mech-generic']
        stress_name = this_parsing_dict['stress_columnn_name']
        strain_name = this_parsing_dict['strain_column_name']
        stress_name_exclusions = this_parsing_dict['stress_exclusion_terms']
        strain_name_exclusions = this_parsing_dict['strain_exclusion_terms']

        #Handle XLSX
        if mech_data_filepath.lower().endswith('.xlsx'):
            #Open the CSV and find the column name row
            excel_df = pd.read_excel(mech_data_filepath, sheet_name= good_sheetname)
            rows = excel_df.values.tolist()
            for row_num, row in enumerate(rows):
                if row_num == start_row_idx:
                    header_row = row
            
            #Walk through header row and find which column index has the stress and strain columns
            for col_idx, cell_text in enumerate(header_row):
                if stress_name in cell_text:
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in stress_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        stress_col_idx = col_idx
                        stress_col_name = cell_text
                if strain_name in cell_text:
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in strain_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        strain_col_idx = col_idx
                        strain_col_name = cell_text

            #Finally, open the data
            data_df = pd.read_csv(mech_data_filepath, usecols=[stress_col_idx, strain_col_idx], skiprows=(start_row_idx) )

            #Make sure the strain is right-side up (i.e. positive values only), ignore lead-in if not at 0, and set minimum at 0
            try:
                strain_series = data_df[strain_col_name].copy()
            except KeyError:
                #if keyerror, parsing has failed and take a look at the data_df to see what's going on
                print(data_df.head(10))
            
              # get info about strain
            pos_strain_sum = len(strain_series[strain_series>0])
            neg_strain_sum = len(strain_series[strain_series<0])
            this_min_strain = strain_series.min()
            
              # flip strain if negative
            if neg_strain_sum > pos_strain_sum:
                # adj_strain_series = strain_series -this_min_strain
                # adj_strain_series = abs(this_min_strain)-adj_strain_series
                adj_strain_series = strain_series * -1
                strain_series = adj_strain_series  #reset series
            
              # check if a lead-in is artificially skewing 0
            data_start_idx = round(len(strain_series)*0.25)  #ignore first part
            new_min = strain_series[data_start_idx::].min()  #find real '0' without lead-in garbage
            strain_series = strain_series + abs(new_min)  #reset series minimum to a more useful 0

              # re-assign strain series
            data_df[strain_col_name] = strain_series

            #Adjust for any 0-stress lead-in
  

            #Assign everything to the output dict
            file_output_dict['dataframe'] = data_df
            file_output_dict['parse_status']['file_pandas_readable'] = True
            file_output_dict['parse_status']['namerow_good'] = True
            file_output_dict['parse_status']['namerow_idx'] = start_row_idx
        
        #Provide cross-compatibility just in case
        elif mech_data_filepath.lower().endswith('.csv'):
            file_output_dict = parse_mech_data_fromcsv(mech_data_filepath,
                                                       data_dict = data_dict)

    return file_output_dict


def get_mech_data_for_name(filename_dict,
                          target_directory = None,
                          do_full_walk = False):
    '''
    Description: Take a 'filename_dict'-like object and return the available data
    '''

    pass


def try_and_guess_units(column_name):
    '''
    Description: Take a string used as a column name and try to extract any stress, strain, or related mechanical units.
    '''

    #Initialize variables
    unit_guess = 'UNK_Units'
    unit_keys = list(default_units_dict.keys())
    good_matches = []
    return_options = []

    #Quick pass through to try and find exact matches
    for unit_string in unit_keys:
        if unit_string in column_name:
            unit_dict = default_units_dict[unit_string]
            if (unit_dict['alternates'] == []):
                if unit_string not in good_matches:
                    good_matches.append(unit_string)
                if unit_dict['exact_match'] not in return_options:
                    return_options.append(unit_dict['exact_match'])
            else:
                for alternate_unit in unit_dict['alternates']:
                    if unit_string not in good_matches:
                        good_matches.append(unit_string)

    #Cull options and select a final guess
    #TODO: add this functionality   

    return unit_guess


def check_logbook_for_printname(printname, logbook):
    '''
    Description:
        Take filename or explicit part print name and return a dictionary of matching printnames from a logbook.
    '''
    #Initialize variables
    matching_rows_dict = {}

    #Parse 'logbook' to make sure it's a DataFrame
    if type(logbook) == str:
        if os.path.isfile(logbook) and logbook.lower().endswith('.csv'):
            logbook = pd.read_csv(logbook)
        else:
            print("Passed string is not a filepath or filepath-like object. Please select an appropriate directory containing a logbook.")
            root = Tk()
            directory = filedialog.askdirectory(title="Select a directory where a logbook copy exists.")
            root.destroy()
            logbook_filepath = get_latest_logbook(directory)
            if logbook_filepath.lower().endswith():
                logbook = pd.read_csv(logbook)

    if isinstance(logbook, pd.DataFrame):
        pass  #hurray

    #Compare 'printname' to printnames in the logbook
      # clean the 'printname' if appropriate


      # get column with printnames from the logbook
    columns = list(logbook.columns)
    possible_printname_columns = []
    for column_name in columns:
        if ('Name' in column_name):
            possible_printname_columns.append(column_name)
    if len(possible_printname_columns) == 1:
        printname_column = possible_printname_columns[0]
    else:
        #TODO: parse different 'name' options to find the printname column
        pass

      # match printname to rows
    matching_row_indices = []
        # attempt exact match
    printname_match_mask = logbook['printname_column'].str.contains(printname, na=False)  #attempt to handle 'nan' exceptions
    matching_rows_dict = logbook[printname_match_mask]

    return matching_rows_dict


def get_latest_logbook(directory):
    '''
    Directory:
        Run through a directory and take a guess at which logbook to use. 
        Format for latest version is "Logbook_AutomatedAnalysisCopy_YYYY-MM-DD". Should be CSV.
    '''
    
    #Initialize variables
    logbook_filepath = None
    potential_logbook_dicts = []
    
    #Walk the directory and try to find a logbook
    for root, dirs, files in os.walk(directory, topdown = True):
        for file in files:

            this_dict = {
                'file': '',
                'filepath': '',
                'immediate_parent_dir': '',
                'filetype': ''
                }

            if ('logbook' in file.lower()) and ('automatedanalysis' in file.lower()):
                logbook_filepath = os.path.join(root, file)

            elif ('logbook' in file.lower()):
                #Instantialize variables
                this_dict['file'] = file
                this_dict['filepath'] = os.path.join(root, file)
                this_dict['immediate_parent_dir'] = os.path.basename(root)

                if file.lower().endswith('.csv'):
                    this_dict['filetype'] = '.csv'
                elif file.lower().endswith('.xlsx'):
                    this_dict['filetype'] = '.xlsx'

                potential_logbook_dicts.update({file: this_dict})

            if not logbook_filepath:
                #TODO: look through potential results and select the 'latest' option
                pass

    return logbook_filepath


def split_filename_for_printname_guessing(filename,
                                          printname_examples = None):
    """
    Description:
        Lorem

    INPUT:
        'filename'          str filepath
    ACTION:
        -lorem
    OUTPUT:
        'printname_dict'    lorem
        
    """

    #Initialize variables
    printname_dict = {
        'full_input_name': '',
        'best_guess_printname': '',
        'iteration_marker': '',

        }


