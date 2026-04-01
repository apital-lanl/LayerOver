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

Description: 

"""

#Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tkinter import filedialog, Tk


  # import other LayerOver modules or module parts
from LayerOver.PSPP.DIWStructure import blank_diw_logbook_row_dict
from LayerOver.PSPP.DIWStructure import standard_logbook_columnnames
from LayerOver.PSPP.DIWStructure import parse_structure

#Define variables
#Define hard-coded thresholds and setting values
default_stress_threshold = 0.2   #in kPa; for silicone elastomers, but should be relatively general
strain_minimum_mask_threshold = -0.1  #minimum strain to accept (<0 to allow for noise at 0 strain)
  # default name for the digital logbook sheet with all the actual logbook data
default_digital_logbok_sheetname = 'Digital Logbook'   #Appropriate sheet as of 2026-01-21

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
    'structure_dict': {}
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

def get_latest_logbook(directory):
    '''
    Directory:
        Run through a directory and take a guess at which logbook to use. 
        Format for latest version is "Logbook_AutomatedAnalysisCopy_YYYY-MM-DD". Should be CSV.
    '''
    
    #Initialize variables
    logbook_filepath = None
    potential_logbook_dicts = {}
    automatedanalysis_versions = []
    
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
                automatedanalysis_versions.append(os.path.join(root, file))

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
    
    best_latest_logbook_guess = None
    if len(automatedanalysis_versions) ==1:
        #Congrats! You found the only appropriate version
        pass
    #If no single 'automatedanalysis' logbook is found, look through the options
    #TODO: finish this
    elif len(automatedanalysis_versions) >1:
        for potential_logbook_filename in list(potential_logbook_dicts.keys()):
            underscore_split = potential_logbook_filename.split('_')
        if len(underscore_split) >1:
            pass
    
    #If there's still no filepath, loop through all found logbooks and try to ID the most appropriate one (latest version)
    elif not best_latest_logbook_guess:
        #TODO: look through potential results and select the 'latest' option
        pass

    return logbook_filepath


def open_logbook(logbook_filepath,
                 target_excel_sheetname = None):
    """
    Description:
        Open a logbook copy and return a DataFrame 

    INPUT:
        'logbook_filepath'  str (filepath); 
    ACTION:
        -lorem
    OUTPUT:
        'logbook_df'        pd.DataFrame; 'cleaned' logbook with standardized column data and column names
    """

    #Initialize variables
    logbook_df = pd.DataFrame([])
    if not target_excel_sheetname:
        target_sheetname =  default_digital_logbok_sheetname

    #Open the logbook
    if logbook_filepath.lower().endswith('.csv'):
        logbook_df = pd.read_csv(logbook_filepath, encoding_errors= 'ignore')
        header_row = list(logbook_df.columns)
    if logbook_filepath.lower().endswith('.xlsx'):
        try:
            logbook_df = pd.read_excel(logbook_filepath, sheet_name= target_sheetname)
        except:
            #TODO: add functionality to look for appropriate sheet if 'target_sheetname' fails
            pass
        header_row = list(logbook_df.columns)
    
    #Clean the dataframe
      # get the column names for columns to be cleaned
    for column_name in header_row:
        #Set some default column names

        if 'strand diameter' in column_name.lower():
            diameter_column_name = column_name
        if 'pitch' in column_name.lower():
            if 'list' in column_name.lower():
                pitch_list_column_name = column_name
            else:
                pitch_column_name = column_name

      # split out skin vs. layer nozzle size if different
      # Logbook column name as of 2026-01-22 = "Strand Diameter, nominal (skin/heli)"
    logbook_df['Strand Diameter, Skin'] = [entry[0] if (len(entry)>1) else entry[0] for entry in logbook_df[diameter_column_name].apply(lambda s: str(s).split(r'/'))]
    logbook_df['Strand Diameter, Layer'] = [entry[1] if (len(entry)>1) else entry[0] for entry in logbook_df[diameter_column_name].apply(lambda s: str(s).split(r'/'))]
    logbook_df['Strand Diameter, Skin'] = logbook_df['Strand Diameter, Skin'].astype(float)
    logbook_df['Strand Diameter, Layer'] = logbook_df['Strand Diameter, Layer'].astype(float)

      # if pitch column as an 'x', replace cell value with layer nozzle size times 'x' ammount (i.e. '1.25x' becomes "1.25 * layer_nozzle_size")
    logbook_df[pitch_column_name] = logbook_df[pitch_column_name].astype(str)
    x_mask = logbook_df[pitch_column_name].str.contains('x', case=False, na=False)
    logbook_df.loc[x_mask, pitch_column_name] = (
        logbook_df.loc[x_mask, pitch_column_name].str.replace('x', '', case=False).astype(float) * 
        logbook_df.loc[x_mask, 'Strand Diameter, Layer']
        )

    logbook_df.loc[x_mask, pitch_list_column_name] = logbook_df.loc[x_mask, pitch_column_name]

    return logbook_df


def check_logbook_for_printname(printname, logbook,
                                row_dict_template = None):
    '''
    Description:
        Take filename or explicit part print name and return a dictionary of matching printnames from a logbook.
    '''
    #Initialize variables
    matching_rows_dict = {}
    if not row_dict_template:
        blank_row_dict = blank_diw_logbook_row_dict
    else:
        blank_row_dict = row_dict_template

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


#######################################################################################################################
#####  Short functions for .apply() on Pandas DataFrame logbook objects  ##############################################
#######################################################################################################################


def add_layernumber_column(logbook_df):
    '''
    Description:
        Take generic 'structure' column from logbook and derive data columns like 'number of layers'.

    INPUT:
        'logbook_df'      pandas.DataFrame; logbook DataFrame object 
    ACTION:
        -lorem
    OUTPUT:
        'logbook_df'      pandas.DataFrame; updated logbook DataFrame
    '''
    def layer_number(row):
        layer_list = parse_structure(row['Structure'])
        layer_length = len(layer_list)
        return layer_length

    logbook_df['Number of Layers'] = logbook_df.apply(layer_number, axis=1)

    return logbook_df

