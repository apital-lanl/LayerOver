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

Description: Module for pulling and analyzing mechanical data (stress/strain exclusively at the moment) from Excel and CSV files.

TODO:
    -'get_latest_logbook'
        -add parsing to try and extract date
        -pull filemetada to get creation and last modified dates
    -'open_logbook'
        -Check 'Pitch' column for 'x' (i.e. "1.25x") and replace with actual value
        -Add "Pitch as strand fraction" column and populate

"""

#Import libraries
import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
from tkinter import filedialog, Tk
import traceback
from scipy.signal import find_peaks

  # import other LayerOver modules or module parts
from LayerOver.PSPP.DIWStructure import blank_diw_logbook_row_dict
from LayerOver.PSPP.DIWStructure import standard_logbook_columnnames
from LayerOver.PSPP.DIWLogbook import open_logbook
from LayerOver.PSPP.DIWLogbook import get_latest_logbook

#Define variables
#Define hard-coded thresholds and setting values
default_stress_threshold = 0.4   #in kPa; for silicone elastomers, but should be relatively general
strain_minimum_mask_threshold = -0.1  #minimum strain to accept (<0 to allow for noise at 0 strain)
  # default name for the digital logbook sheet with all the actual logbook data
default_digital_logbok_sheetname = 'Digital Logbook'   #Appropriate sheet as of 2026-01-21
# assymptotic_stress_threshold = 1000   # in kPa; high number for 'max lock-up stress threshold'; quadratic stress increase is calculated as mostly-linear, so a huge number is needed and isn't physical
assymptotic_stress_max_multiplier = 1.5  # assymptote is drawn and 'multiplier' * 'max_stress'; i.e. at 1.5, a 500 kPa peak would be assymptoted at 750

#Example column names for various types of report from mechanical testing instruments
# used to guess which 1) type of mech data is being parsed, 2) which row in the spreadsheet contains the column names, 
# and 3) which columns to pull for 'stress' and 'strain' data.
example_namerows = {
    'mech-generic': ['Load', 'Stress', 'Strain'],
    'mech-type1': ['Index','Time (sec)', 'Load (kN)', 'Crosshead (mm)', 'PrimaryExtension (mm)', 'Stress (kPa)', 'Strain (mm/mm)', 'Specimen Height (in)'],
    'mech-type2': ['Index','Time (sec)', 'Load (N)', 'Stress (kPa)', 'Gap (mm)', 'GapEx1 (mm)', 'Extension (mm)', 'Extension Ex1 (mm)', 'Strain (mm/mm)', 'Strain Ex1 (mm/mm)'],
    'mech-type3': ['Time (sec)', 'Load (N)', 'Stress (kPa)', 'Uncorrected Gap (mm)', 'Corrected Gap (mm)', 'Strain (mm/mm)']
    }


#Details for parsing files of the types in 'example_namerows'
#'__exclusion_terms' are lists of terms that disqualify a name. I.E. 'Strain Ex1 (mm/mm)' is a valid 'strain' name for 'mech-type2', but we only want the 'Strain (mm/mm)' column and 
#       and should ignore any 'strain' column with ' Ex1' also in the name.
namerow_example_parsing_dict = {
    'mech-generic': {
        'stress_columnn_name': 'Stress',
        'strain_column_name': 'Strain',
        'extension_column_name': 'Extension',
        'stress_exclusion_terms': ['peak'],
        'strain_exclusion_terms': [' Ex1 '],
        'extension_exclusion_terms': []
        },
    'mech-type1': {
        'stress_columnn_name': 'Stress (kPa)',
        'strain_column_name': r'Strain (mm/mm)',
        'extension_column_name': 'PrimaryExtension (mm)',
        'stress_exclusion_terms': ['peak'],
        'strain_exclusion_terms': [],
        'extension_exclusion_terms': []
        },
    'mech-type2': {
        'stress_columnn_name': 'Stress (kPa)',
        'strain_column_name': r'Strain (mm/mm)',
        'extension_column_name': 'Extension (mm)',
        'stress_exclusion_terms': ['peak'],
        'strain_exclusion_terms': [],
        'extension_exclusion_terms': ['ex1']
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

#Stress at the following strain values will be reported in Mecahnical Summary folder and CSV
#Any strain values in 0.01 increments can be called
strain_reporting_points = [
    0.01,
    0.02,
    0.03,
    0.04,
    0.05,
    0.06,
    0.07,
    0.08,
    0.09,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.85,
    0.9
    ]


######################################################################################################################################
##### Generic Utilities  #############################################################################################################
######################################################################################################################################

def process_directory_for_mech_files(directory=None, 
                                     show_each_file_results = False,
                                     save_last_replicate_graph = False):
    
    #Initialize variables
    mech_data_dict = {}
      # if no directory is added as input, select one
    if not directory:
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()
    #Try and find the 'latest' logbook
      # brute-force a logbook search for a couple of levels of parent directory if none is found in analysis directory
    logbook_filepath = get_latest_logbook(directory)   #Returns a str if successful and NONE if not successful
      # immediate parent directory
    if not logbook_filepath:
        parent_directory = os.path.basename(os.path.dirname(directory))
        logbook_filepath = get_latest_logbook(parent_directory)
      # grand-parent directory
    if not logbook_filepath:
        parent_directory = os.path.basename(os.path.dirname(parent_directory))
        logbook_filepath = get_latest_logbook(parent_directory)
      # open the logbook as a pandas DataFrame
    if logbook_filepath:
        logbook_df = open_logbook(logbook_filepath,
                                  target_excel_sheetname = None)
        logbook_df_column_names = list(logbook_df.columns)
    else:
        #TODO: need to add final recourse here if no logbook has been found
        pass

    #Get all CSV and XLSX files in the directory
    mech_filedata_dict = walk_directory_for_mech_files(directory)
    mech_file_keys = list(mech_filedata_dict.keys())

    #Open each file and process
    for file_idx, filename_key in enumerate(mech_file_keys):

        #Initialize this file's return dict
        this_file_dict = {
            'filename': None,
            'filepath': None,
            'filetype': None,
            'immediate_parent_directory': None,
            'parse_type': None,
            'clean_filename': None,
            }

        #Pull file metadata
        filedata_dict = mech_filedata_dict[filename_key]
        data_filepath = filedata_dict['filepath']
          # assign values to output dict
        this_file_dict['filetype'] = filedata_dict['filetype']
        this_file_dict['filename'] = filedata_dict['filename']
        this_file_dict['filepath'] = filedata_dict['filepath']
        this_file_dict['immediate_parent_directory'] = filedata_dict['immediate_parent_directory']
        this_file_dict['parse_type'] = filedata_dict['parse_type']
        this_file_dict['clean_filename'] = filedata_dict['clean_filename']

        print()
        print("#"*50)
        print(f"Parsing {os.path.basename(this_file_dict['filename'])}")
        print(f"\t(file {file_idx+1} of {len(mech_file_keys)})")

        #Take mechanical data and populate a return dictionary from various processing steps
        #  NOTE: these are the steps that sometimes fail, so they're wrapped into a Try loop
        #  TODO: add error reporting functionality 
        try:
            #Process the actual file
            process_dict = process_mechanical_file_against_logbook(data_filepath, logbook_df,
                                                                    alternate_save_directory = directory,
                                                                    show_each_file_results = show_each_file_results,
                                                                    save_last_replicate_graph = save_last_replicate_graph)

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
            this_file_dict['successful_columname_parse'] = process_dict['data_namerow_dict']['successful_parse']
            this_file_dict['file_namerow_index'] = process_dict['data_namerow_dict']['likely_name_row']
            this_file_dict['pandas_readable'] = process_dict['pandas_readable']
            this_file_dict['strain_assymptote'] = process_dict['last_strain_lockup']
            this_file_dict['extension_assymptote'] = process_dict['last_extension_lockup']
            strain_reporting_keys = list(process_dict['stress_dict'].keys())
            for strain_value in strain_reporting_keys:
                this_name = f"stress_at_{strain_value}_strain"
                try:
                    this_file_dict[this_name] = process_dict['stress_dict'][strain_value]
                except:
                    this_file_dict[this_name] = 'nan'
            # this_file_dict[''] = process_dict['']
            # this_file_dict['raw_data'] = data_df  #not sure we actually want to return this, but we could
            
            #Clean the filename and check against logbook
            printname_guess_dict =   split_filename_for_printname_guessing(this_file_dict['filename'], logbook_df)
            matches = printname_guess_dict['all_matches']
            scores = printname_guess_dict['all_match_scores']
            best_guess = printname_guess_dict['best_guess_printname']
            best_iterator_guess = printname_guess_dict['iteration_marker']
            parse_check = printname_guess_dict['printname_parse_bool']
            log_book_entry = printname_guess_dict['logbook_entry']

            #Add results to the return row for this file
            this_file_dict['logbook_entry_found'] = parse_check
            if parse_check:
                this_file_dict['mechanical_data_iteration'] = best_iterator_guess
                for logbook_column in logbook_df_column_names:
                    #Flag that logbooks might be altered by user (i.e. u instead of micron symbol)
                    try:
                        these_logbook_entries = log_book_entry[logbook_column]
                        #Rare, but if a print fails there may be more than one entry for the same print name
                        if these_logbook_entries.shape[0] == 1:
                            this_file_dict[logbook_column] = these_logbook_entries.values[0]
                          # if there's more than one entry for a printname, assume the last entry is the 'good' one
                        else:
                            #TODO: add functionality to actually parse the row for 'failure' and find proper print if it exists
                            this_file_dict[logbook_column] = these_logbook_entries.values[-1] 
                    except:
                        this_file_dict[logbook_column] = ''

            #Write the results to the global directory return dictionary    
            mech_data_dict[filename_key] = this_file_dict

              # report success
            print()
            print("Entry added.")
            print("_"*50)
            print()

        except Exception as le:
            #Write the results to the global directory return dictionary    
            mech_data_dict[filename_key] = this_file_dict

            #'le' stands for 'loop error'
            print(f"Exception thrown: {le}")
            print(traceback.format_exc())
            print()
            print("Failure to parse; moving on to next file.")
            print("_"*50)
            print()
            
    #Make the final dataframe
    #  Add each key:value as a row to the dataframe; key is ignored
    meta_data_df = pd.DataFrame.from_dict(mech_data_dict, orient= "index")

    #Save the final dictionary
    save_directory = os.path.join(directory, "Extracted Mechanical Summary Data")
    csv_save_name = os.path.join(save_directory, "Summary of Directory Mechanical Data.csv")
    meta_data_df.to_csv(csv_save_name, index = False)

    return meta_data_df


def process_mechanical_file_against_logbook(data_filepath, logbook_df,
                                            alternate_save_directory = None,
                                            show_each_file_results = False,
                                            show_displacement_peaks = False,
                                            save_last_replicate_graph = False,
                                            displacement_column = 'extension'):
    """
    Description:
        Open a mechanical data file (CSV or XLSX) and attempt parse data relative to a logbook 
    
    INPUTS:
        'data_filepath'             str; filepath to CSV or XLSX mechanical data.
        'logbook_df'                pd.DataFrame; logbook-like-DataFrame with information for printnames that the mechanical data file can be matched against.
        (OPTIONAL)
        'alternate_save_directory'  str; directory path for a global save locaiton; i.e. if processing a whole directory, that directory path is passed
                                    here to allow for ALL the directory's file results to be saved in one location.
                                    NOTE: ONLY PASS A PRIMARY DIRECTORY. If  directory location is passed, new folders are created in that location. 
        'show_each_file_results'
        'save_last_replicate_graph' lorem;
        'displacement_column'       str; 'strain'  'extension'
            'strain'                use strain as the displacement column
            'extension'             use raw data as the displacement column; i.e. raw extension 0-indexed to first value with real stress
        'displacement_zero_offset'  int; number of indices back from 'zero' to include for strain and extension

    ACTION:
        -lorem
    OUTPUTS:
        'lorem'     lorem
    """
    
    #Initialize and condition variables
    this_file_dict = {
        'filename': None,
        'filepath': None,
        'filetype': None,
        'immediate_parent_directory': None,
        'parse_type': None,
        'clean_filename': None,
        }
    logbook_df_column_names = list(logbook_df.columns)
      # if filepath is bad, return trivial file dictionary
    if not os.path.isfile(data_filepath):
        print(f"Bad filepath; {data_filepath}")
        return this_file_dict
    else:
        #Get standard filepath properties
        this_file_dict['filename'] = os.path.basename(data_filepath).replace('.csv','').replace('.xlsx','')
        this_file_dict['filepath'] = data_filepath
        this_file_dict['immediate_parent_directory'] = os.path.basename(os.path.dirname(data_filepath))
        
        #Get mechanical-data-specific filename properties
        parse_dict = parse_mech_data_filename(data_filepath)  #returns dict with 'raw_filename', 'parse_type', 'clean_filename'
        this_file_dict['parse_type'] = parse_dict['parse_type']
        this_file_dict['clean_filename'] = parse_dict['clean_filename']

    #Load and parse the actual raw data from the file
    if data_filepath.lower().endswith('csv'):
        file_output_dict = parse_mech_data_fromcsv(data_filepath)
        data_df = file_output_dict['dataframe']
        this_file_dict['filetype'] = 'csv'
    elif data_filepath.lower().endswith('xlsx'):
        namerow_dict = check_for_namerow(data_filepath, show_peaks = False)
        good_sheet_check = bool(namerow_dict['successful_parse'] and (namerow_dict['data_sheetname'] != ''))
        if good_sheet_check:
            file_output_dict = parse_mech_data_fromxlsx(data_filepath, data_dict=namerow_dict)
            data_df = file_output_dict['dataframe']
            this_file_dict['filetype'] = 'xlsx'
        else:
            print("Failure to find good sheet in Excel file; check parsing or file contents.")
            data_df = pd.DataFrame({"Failure":[]})

    #Validate namerow location for column names
    namerow_dict = check_for_namerow(data_filepath, 
                                     show_peaks = show_each_file_results)
        # assign values to output dict
    this_file_dict['data_namerow_dict'] = namerow_dict
        
    #Assign data and parse status to output dict
    this_file_dict['pandas_readable'] = file_output_dict['parse_status']['file_pandas_readable']
    this_file_dict['raw_data'] = data_df  #not sure we actually want to return this, but we could

    #Take mechanical data and populate a return dictionary from various processing steps
    #  NOTE: these are the steps that sometimes fail, so they're wrapped into a Try loop
    #  TODO: add error reporting functionality 
    try:
        #Assign data and find correct columns
        data_df = file_output_dict['dataframe']
        column_names = list(data_df.columns)
        for column in column_names:
            if 'stress' in column.lower():
                stress_col_name = column
            if 'strain' in column.lower():
                strain_col_name = column
            if 'extension' in column.lower():
                extension_col_name = column

        print()
        print(data_df.head(7))

        #Plot all the curves
        if show_each_file_results:
            plt.figure(figsize = (10,10))
            if displacement_column == 'strain':
                plt.scatter(data_df[strain_col_name], data_df[stress_col_name])
                plt.xlabel("Strain (mm/mm)")
            elif displacement_column == 'extension':
                plt.scatter(data_df[extension_col_name], data_df[stress_col_name])
                plt.xlabel("Extension (mm)")
            plt.title(f"All mech data in {os.path.basename(data_filepath)}")
            plt.ylabel("Stress")
            plt.show()

        replicate_dict = pull_mechanical_replicates(data_df, data_dict = None, 
                                                       show_peaks = show_displacement_peaks,
                                                       stress_threshold = None,
                                                       displacement_zero_offset = 1,
                                                       strain_min_thresh = None,
                                                       report_nonnegative_strain = True,
                                                       displacement_column = displacement_column,
                                                       points_n_for_assymptote = 5)
                                   
            # pull keys from returned dict
        replicate_parse_success =  replicate_dict['replicate_parse_success']
        number_of_replicates = replicate_dict['number_of_replicates']
        last_strain_peak_index = replicate_dict['last_strain_peak_index']
        last_strain_valley_index = replicate_dict['last_strain_valley_index']
        peak_to_valley_index_diff = replicate_dict['peak_to_valley_index_diff']
        #The following are passed as DataFrames because I was too lazy to change it
        #Convert to values and just pull the first one; otherwise a series is returned with N identical values
        last_strain_lockup = replicate_dict['replicate_data'][number_of_replicates]['strain_lockup'].values[0]
        this_file_dict['last_strain_lockup']= last_strain_lockup
        last_extension_lockup = replicate_dict['replicate_data'][number_of_replicates]['extension_lockup'].values[0]
        this_file_dict['last_extension_lockup']= last_extension_lockup
        this_file_dict['stress_dict']= replicate_dict['last_replicate_stress_dict']

            # each value in "replicate_data_dict" is a pandas.DataFrame (hopefully)
            # replicate numbering starts at 1
        replicate_data_dict = replicate_dict['replicate_data']
        this_file_dict['replicate_data'] = replicate_dict['replicate_data']   #not sure we actually want to return this, but we could

        #Try and pull, plot the last replicate if desired
        try:
            last_df = replicate_data_dict[number_of_replicates]

            print()
            print(" "*5, "#"*15)
            print()
            if replicate_parse_success:
                print(f"Number of replicates parsed: {number_of_replicates}")
            else:
                print("Failure to resolve individual mechanical data replicates.")
            print()
        except KeyError:
            print()
            print("Failure to resolve individual mechanical data replicates")
            print("     (no replicate mechanical data found for last index)")
        try:
            if show_each_file_results:
                if displacement_column == 'strain':
                    plt.figure(figsize = (10,10))
                    plt.scatter(last_df['strain_data_loading'], last_df['stress_data_loading'], color = 'r')
                    plt.scatter(last_df['strain_data_unloading'], last_df['stress_data_unloading'], color = 'g')
                    plt.plot([last_strain_lockup, last_strain_lockup], [0,last_df['stress_data_loading'].max()], color = 'lime', linewidth = 4)
                    plt.title(f"Final replicate cycle for {os.path.basename(data_filepath)}")
                    plt.xlabel("Strain")
                    plt.ylabel("Stress")
                    plt.legend([f"Replicate {number_of_replicates} Loading curve", f"Replicate {number_of_replicates} Unloading curve", "Lockup assymp. guess"])
                elif displacement_column == 'extension':
                    plt.figure(figsize = (10,10))
                    plt.scatter(last_df['extension_data_loading'], last_df['stress_data_loading'], color = 'r')
                    plt.scatter(last_df['extension_data_unloading'], last_df['stress_data_unloading'], color = 'g')
                    plt.plot([last_extension_lockup, last_extension_lockup], [0,last_df['stress_data_loading'].max()], color = 'lime', linewidth = 4)
                    plt.title(f"Final replicate cycle for {os.path.basename(data_filepath)}")
                    plt.xlabel("Extension")
                    plt.ylabel("Stress")
                    plt.legend([f"Replicate {number_of_replicates} Loading curve", f"Replicate {number_of_replicates} Unloading curve", "Lockup assymp. guess"])
                if save_last_replicate_graph:
                    if alternate_save_directory:
                        plot_save_directory = os.path.join(alternate_save_directory, "Extracted Mechanical Summary Data", "Last Stress-Strain Replicate Graphs")
                        if not os.path.isdir(plot_save_directory):
                            os.mkdir(os.path.join(alternate_save_directory, "Extracted Mechanical Summary Data"))
                            os.mkdir(os.path.join(alternate_save_directory, "Extracted Mechanical Summary Data", "Last Stress-Strain Replicate Graphs"))
                    else:
                        plot_save_directory = os.path.dirname(data_filepath)
                    plot_save_name = os.path.join(plot_save_directory, f"LastMechReplicateGraph_{this_file_dict['filename']}.png")
                    plt.savefig(plot_save_name, dpi=300)
                plt.show()
            elif save_last_replicate_graph:
                if displacement_column == 'strain':
                    plt.figure(figsize = (10,10))
                    plt.scatter(last_df['strain_data_loading'], last_df['stress_data_loading'], color = 'r')
                    plt.scatter(last_df['strain_data_unloading'], last_df['stress_data_unloading'], color = 'g')
                    plt.plot([last_strain_lockup, last_strain_lockup], [0,last_df['stress_data_loading'].max()], color = 'lime', linewidth = 4)
                    plt.title(f"Final replicate cycle for {os.path.basename(data_filepath)}")
                    plt.xlabel("Strain")
                    plt.ylabel("Stress")
                    plt.legend([f"Replicate {number_of_replicates} Loading curve", f"Replicate {number_of_replicates} Unloading curve", "Lockup assymp. guess"])
                elif displacement_column == 'extension':
                    plt.figure(figsize = (10,10))
                    plt.scatter(last_df['extension_data_loading'], last_df['stress_data_loading'], color = 'r')
                    plt.scatter(last_df['extension_data_unloading'], last_df['stress_data_unloading'], color = 'g')
                    plt.plot([last_extension_lockup, last_extension_lockup], [0,last_df['stress_data_loading'].max()], color = 'lime', linewidth = 4)
                    plt.title(f"Final replicate cycle for {os.path.basename(data_filepath)}")
                    plt.xlabel("Extension")
                    plt.ylabel("Stress")
                    plt.legend([f"Replicate {number_of_replicates} Loading curve", f"Replicate {number_of_replicates} Unloading curve", "Lockup assymp. guess"])
                if alternate_save_directory:
                    plot_save_directory = os.path.join(alternate_save_directory, "Extracted Mechanical Summary Data", "Last Stress-Strain Replicate Graphs")
                    if not os.path.isdir(plot_save_directory):
                        os.mkdir(os.path.join(alternate_save_directory, "Extracted Mechanical Summary Data"))
                        os.mkdir(os.path.join(alternate_save_directory, "Extracted Mechanical Summary Data", "Last Stress-Strain Replicate Graphs"))
                else:
                    plot_save_directory = os.path.dirname(data_filepath)
                plot_save_name = os.path.join(plot_save_directory, f"LastMechReplicateGraph_{this_file_dict['filename']}.png")
                plt.savefig(plot_save_name, dpi=300)
                plt.close()

        except:
            print()
            print("Final replicate plot failure (hopefully for obvious reasons)")
            print(traceback.format_exc())
            print()
                
        #Clean the filename and check against logbook
        printname_guess_dict =   split_filename_for_printname_guessing(this_file_dict['filename'], logbook_df)
        matches = printname_guess_dict['all_matches']
        scores = printname_guess_dict['all_match_scores']
        best_guess = printname_guess_dict['best_guess_printname']
        best_iterator_guess = printname_guess_dict['iteration_marker']
        parse_check = printname_guess_dict['printname_parse_bool']
        log_book_entry = printname_guess_dict['logbook_entry']

        #Add results to the return row for this file
        this_file_dict['logbook_entry_found'] = parse_check
        if parse_check:
            this_file_dict['mechanical_data_iteration'] = best_iterator_guess
            for logbook_column in logbook_df_column_names:
                #Flag that logbooks might be altered by user (i.e. u instead of micron symbol)
                try:
                    these_logbook_entries = log_book_entry[logbook_column]
                    #Rare, but if a print fails there may be more than one entry for the same print name
                    if these_logbook_entries.shape[0] == 1:
                        this_file_dict[logbook_column] = these_logbook_entries.values[0]
                      # if there's more than one entry for a printname, assume the last entry is the 'good' one
                    else:
                        #TODO: add functionality to actually parse the row for 'failure' and find proper print if it exists
                        this_file_dict[logbook_column] = these_logbook_entries.values[-1] 
                except:
                    this_file_dict[logbook_column] = ''
        else:
            this_file_dict['mechanical_data_iteration'] = best_iterator_guess
            for logbook_column in logbook_df_column_names:
                this_file_dict[logbook_column] = ''

        #Write the results to the global directory return dictionary    
        return this_file_dict

    except Exception as le:
        #'le' stands for 'loop error'
        print(f"Exception thrown: {le}")
        print(traceback.format_exc())
        print()
        print(f"Failure to parse {data_filepath}")
        print("_"*50)
        print()

        return this_file_dict


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
    Description: 
        Generic function to open a spreadsheet (CSV or XLSX), find likely namerow based on an example or largest text-containing row, 
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
            namerow_dict['data_sheetname'] = sheet_name
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
                                   show_peaks = False,
                                   stress_threshold = None,
                                   displacement_zero_offset = 1,
                                   strain_min_thresh = None,
                                   report_nonnegative_strain = True,
                                   displacement_column = 'extension',
                                   points_n_for_assymptote = 5):
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
    #Condition variables
    displacement_zero_offset = int(round(displacement_zero_offset))  #make sure zero_offset is an integer
    strain_zero_offset = displacement_zero_offset
    extension_zero_offset = displacement_zero_offset

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
        'all_extension_data': []
        }

    data_dict = {
        'mech_units_info': {
            'testing_index': 'UNK',
            'stress_units': 'UNK',
            'strain_units': 'UNK',
            'extension_units': 'UNK'
            },
        }
    columns = list(data_df.columns)
    stress_flag = False
    strain_flag = False
    extension_flag = False

    #Check for 'data_df' data
    #TODO: add DataFrame type check

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

        if ('extension' in column_name.lower()) and not extension_flag:
            mech_dict['all_extension_data'] = data_df[column_name]
            print()
            print(f"Added {column_name} as extension data")
            extension_flag = True
            data_dict['mech_units_info']['extension_units'] = unit_guess
            extension_units = unit_guess
            extension_col_name = f"Extension ({extension_units})"
        elif ('extension' in column_name.lower()):
            print
            print(f"Multiple extension columns passed. Ignoring {column_name}")

    #Create mask of high-stress and low-stress regions
    stress_max = mech_dict['all_stress_data'].max()
    high_stress_mask = mech_dict['all_stress_data'].copy() > stress_max * 0.9
    low_stress_mask = mech_dict['all_stress_data'].copy() < 5   #in kPa

    #Zero-correct (eliminate negative values) and max-correct raw data
    #  Should flip the x-axis values (bigger magnitude -> smaller magnitude)
      # raw extension values (typically mm)
    ext_pos_check = bool(mech_dict['all_extension_data'][high_stress_mask][-10:-1].sum() > mech_dict['all_extension_data'][low_stress_mask][-10:-1].sum())
    if not ext_pos_check:   #'False' if low stress extension values are bigger than the peak stress extension values
        initial_max = mech_dict['all_extension_data'].max()
        mech_dict['all_extension_data'] = mech_dict['all_extension_data'] - initial_max
        mech_dict['all_extension_data'] = mech_dict['all_extension_data'].abs()
    #   # normalized strain values (mm/mm)
    strain_pos_check = bool(mech_dict['all_strain_data'][high_stress_mask][-10:-1].sum() > mech_dict['all_strain_data'][low_stress_mask][-10:-1].sum())
    if not strain_pos_check:   #'False' if low stress extension values are bigger than the peak stress extension values
        initial_max = mech_dict['all_strain_data'].max()
        mech_dict['all_strain_data'] = mech_dict['all_strain_data'] - initial_max
        mech_dict['all_strain_data'] = mech_dict['all_strain_data'].abs()

    #Use peak finding to grab last loading/unloading cycl
    raw_df = pd.DataFrame(data= {stress_col_name: mech_dict['all_stress_data'],
                                  strain_col_name: mech_dict['all_strain_data'],
                                  extension_col_name: mech_dict['all_extension_data']})
    
      # get reasonable peak height and distance between peak expectations
    length = raw_df[f'Strain ({strain_units})'].values.shape[0]
    distance_guess = length * 0.15
    height_guess = raw_df[f'Strain ({strain_units})'].values[round(length*0.25)::].max()
      # find all peaks in strain data
      # NOTE: both distance and height are imperfect proxies for wonky strain peaks, but distance tends to perform better in avoiding multiple non-peak values
      #       'distance' in this case is an expectation of how far apart peak values should be, set by 'distance_guess' above
    # peaks, _ = find_peaks(raw_df[f'Strain ({strain_units})'].values, height = height_guess)
    peaks, _ = find_peaks(raw_df[f'Strain ({strain_units})'].values, distance = distance_guess)
    replicate_dict['number_of_replicates'] = len(peaks)
    
    if len(peaks) == 0:
        print("No peaks found in strain data")
    #Use strain data to define replicate regions of the entire raw data Series
    elif len(peaks) == 1:
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
        strain_max = raw_df[strain_col_name].max()
        extension_max = raw_df[extension_col_name].max()

        #Generate the graph
        if displacement_column == 'strain':
            plt.figure(figsize = (10,10))
            plt.scatter(list(range(raw_df[strain_col_name].shape[0])), raw_df[strain_col_name])
            plt.title(f"Strain peaks found")
            #Show peak start/stop locations defining replicates
            for peak in peaks:
                plt.scatter(peak, raw_df[strain_col_name][peak], marker = 'x', s = 200, color = 'gray')
              # start of last loading replicate
            plt.plot([last_valley_guess_index, last_valley_guess_index], [0, strain_max], color='r', linewidth=3, alpha=0.6)
              # start of last unloading replicate
            plt.plot([last_peak_index, last_peak_index], [0, strain_max], color='g', linewidth=3, alpha=0.6)
            plt.xlabel("Index of strain value")
            plt.ylabel("Strain value (mm/mm)")
            plt.show()
        elif displacement_column == 'extension':
            plt.figure(figsize = (10,10))
            plt.scatter(list(range(raw_df[extension_col_name].shape[0])), raw_df[extension_col_name])
            plt.title(f"Extension peaks found")
            #Show peak start/stop locations defining replicates
            for peak in peaks:
                plt.scatter(peak, raw_df[extension_col_name][peak], marker = 'x', s = 200, color = 'gray')
              # start of last loading replicate
            plt.plot([last_valley_guess_index, last_valley_guess_index], [0, extension_max], color='r', linewidth=3, alpha=0.6)
              # start of last unloading replicate
            plt.plot([last_peak_index, last_peak_index], [0, extension_max], color='g', linewidth=3, alpha=0.6)
            plt.xlabel("Index of extension value")
            plt.ylabel("Extension value (mm)")
            plt.show()
    #If replicate parsing not successful, show a graph that might hint at why
    elif show_peaks:
        plt.figure(figsize = (10,10))
        if displacement_column == 'strain':
            plt.scatter(list(range(raw_df[strain_col_name].shape[0])), raw_df[strain_col_name])
            plt.xlabel("Strain (mm/mm)")
        elif displacement_column == 'extension':
            plt.scatter(list(range(raw_df[extension_col_name].shape[0])), raw_df[extension_col_name])
            plt.xlabel("Extension (mm)")
        plt.ylabel("Stress")
        plt.title(f"Failed to parse replicates from this data")
        plt.show()
    
    #Loading cycles start at (peak_idx-valley_idx), Unloading cycles end at (peak_idx+valley_idx)
    for replicate_idx, peak_idx in enumerate(peaks):
        # load_start_idx = peak_idx-cycle_index_diff_guess-strain_zero_offset   #too many issues with this offset
        load_start_idx = peak_idx-cycle_index_diff_guess
        unload_end_idx = peak_idx+cycle_index_diff_guess
        #Make sure indices are in-bounds for data size
        if load_start_idx<0:
            load_start_idx = 0
          # should be the same size, but just in case take the minimum of the two columns dimensions
        smallest_data_index = min(raw_df[strain_col_name].shape[0], raw_df[stress_col_name].shape[0])-1  #subtract 1 for 0-start indexing
        if unload_end_idx > smallest_data_index:  #make sure index is valid
            unload_end_idx = smallest_data_index

        #Find strain minimum where stress actually starts increasing above a threshold
        #We'll only consider the first loading cycle
          ## Use first 'replicate' to avoid lead-in garbage data
        if (replicate_idx == 0):
            #If no setting is passed, use module default
            if not stress_threshold:
                stress_threshold = default_stress_threshold

            #Grap the loading curve for analysis of displacement 0 point relative to when stress begins to show
            strain_loading = raw_df[strain_col_name].iloc[load_start_idx:peak_idx]
            stress_loading = raw_df[stress_col_name].iloc[load_start_idx:peak_idx]
            extension_loading = raw_df[extension_col_name].iloc[load_start_idx:peak_idx]

            #Find stress index when stress crosses threshold value (hard-coded in module defaults; in kPa typically)
            stress_loading_avg = stress_loading.rolling(5, center=True, min_periods = 1).mean()
            valid_stress_mask = stress_loading_avg >= stress_threshold
            first_valid_stress_index = valid_stress_mask.idxmax()-4  #returns first True value; subtract 4 to account for the rolling average window size (5) and center=True
            first_valid_strain_index = first_valid_stress_index-strain_zero_offset   #include a few prior readings to get a clear lead-in to the stress; not a big deal typically
            if first_valid_strain_index < 0:
                first_valid_strain_index = 0
              ## duplicate strain index start for extension; set separately here to allow for specific extension-related compensation in the future
            first_valid_ext_index = first_valid_strain_index

            #Get displacement values at which stress actually starts to rise to normalize displacement values to
              ## get these values before resetting the strain data so we can apply the same offset to all replicates
            first_valid_strain_value = strain_loading.iloc[first_valid_strain_index]
            # first_valid_strain_value = strain_loading.values[first_valid_strain_index]  #alternative
            first_valid_extension = extension_loading.iloc[first_valid_strain_index]
            # first_valid_extension = extension_loading.values[first_valid_ext_index]  #alternative
              ## get the strain value at new predicted '0 strain' value
            strain_zero_correction = first_valid_strain_value
            extension_zero_correction = first_valid_extension

            #Apply 0-displacement corrections to the entire series
            raw_df[strain_col_name] = raw_df[strain_col_name]- strain_zero_correction
            raw_df[extension_col_name] = raw_df[extension_col_name]- extension_zero_correction

        #Grab cycle data
          # pull cyle data
        strain_loading = raw_df[strain_col_name].iloc[load_start_idx:peak_idx]
        stress_loading = raw_df[stress_col_name].iloc[load_start_idx:peak_idx]
        extension_loading = raw_df[extension_col_name].iloc[load_start_idx:peak_idx]
        strain_unloading = raw_df[strain_col_name].iloc[peak_idx:unload_end_idx]
        stress_unloading = raw_df[stress_col_name].iloc[peak_idx:unload_end_idx]
        extension_unloading = raw_df[extension_col_name].iloc[peak_idx:unload_end_idx]
        
        #If flagged, make sure non-negative data is reported
        if report_nonnegative_strain:
            #If no min strain threshold is passed, use module default
            if not strain_min_thresh:
                strain_min_thresh = strain_minimum_mask_threshold
            #Get pandas Series mask for values of strain above threshold
            loading_mask = strain_loading>= 0
            unloading_mask = strain_unloading>= 0
            #Apply non-negative mask to get appropriate values only
            strain_loading = strain_loading[loading_mask]
            stress_loading = stress_loading[loading_mask]
            strain_unloading = strain_unloading[unloading_mask]
            stress_unloading = stress_unloading[unloading_mask]

        #Try and find an assymptote for the loading curve
          # make a stress threshold to assymptote to
        assymptotic_stress_threshold = stress_loading.max()*assymptotic_stress_max_multiplier
          # get derivatives
          # only pull the positive derivates to avoid issues from the way data is collected (probably because of time resolution errors and noise)
        extension_loading_diff = extension_loading.diff()
        extension_loading_diff = extension_loading_diff[extension_loading_diff>0]
        strain_loading_diff = strain_loading.diff()
        strain_loading_diff = strain_loading_diff[strain_loading_diff>0]
        loading_stress_diff = stress_loading.diff()
        loading_stress_diff = loading_stress_diff[loading_stress_diff>0]
          # get second derivatives
        extension_loading_d2 = extension_loading_diff.diff()
        extension_loading_d2 = extension_loading_d2[extension_loading_d2>0]
        strain_loading_d2 = strain_loading_diff.diff()
        strain_loading_d2 = strain_loading_d2[strain_loading_d2>0]
        loading_stress_d2 = loading_stress_diff.diff()
        loading_stress_d2 = loading_stress_d2[loading_stress_d2>0]
        #Prior version
          # pull last 10 values and get a slope for the lock-up acceleration
        if (extension_loading_d2.tail(points_n_for_assymptote).mean()) < 1e-6:    
            extension_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / 1e-6   #avoid division by zero with a very low threshold
        else:
            extension_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / (extension_loading_d2.tail(points_n_for_assymptote).mean())
        if (strain_loading_d2.tail(points_n_for_assymptote).mean()) < 1e-6:    
            strain_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / 1e-6   #avoid division by zero with a very low threshold
        else:
            strain_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / (strain_loading_d2.tail(points_n_for_assymptote).mean())

        # pull last 10 values and get a slope for the lock-up acceleration
        if points_n_for_assymptote < 1e-6:    
            extension_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / 1e-6   #avoid division by zero with a very low threshold
        else:
            extension_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / points_n_for_assymptote
        if points_n_for_assymptote < 1e-6:    
            strain_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / 1e-6   #avoid division by zero with a very low threshold
        else:
            strain_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / points_n_for_assymptote
        stress_accel = (loading_stress_d2.tail(points_n_for_assymptote).mean()) / points_n_for_assymptote
          #strain assymptote
          #NOTE: could solve  analytically, but this is easier and more clear; never more than a few thousand of these to do in total
        last_stress = stress_loading.tail(1).values[0]
        last_strain = strain_loading.tail(1).values[0]
        strain_step = strain_loading_diff.mean()
        strain_velocity = strain_loading_diff.tail(3).mean()
        stress_velocity = loading_stress_diff.tail(3).mean()
        iteration_counter = 0
        strain_delta = 0
        while (last_stress<assymptotic_stress_threshold) and (iteration_counter <100):
            # last_stress = last_stress + stress_velocity*(strain_step) + 0.5 * stress_accel*(strain_step)**2
            last_stress = last_stress + stress_velocity + 0.5 * stress_accel
            strain_delta = strain_delta + strain_step
            iteration_counter += 1
        strain_lockup = last_strain + strain_delta
          #extension assymptote
          #NOTE: could solve  analytically, but this is easier and more clear; never more than a few thousand of these to do in total
        last_stress = stress_loading.tail(1).values[0]
        last_ext = extension_loading.tail(1).values[0]
        ext_step = extension_loading_diff.mean()
        ext_velocity = extension_loading_diff.tail(3).mean()
        iteration_counter = 0
        ext_delta = 0
        while (last_stress<assymptotic_stress_threshold) and (iteration_counter <100):
            # last_stress = last_stress + stress_velocity*(ext_step) + 0.5 * stress_accel*(ext_step)**2
            last_stress = last_stress + stress_velocity + 0.5 * stress_accel
            ext_delta = ext_delta + ext_step
            iteration_counter += 1
        ext_lockup = last_ext + ext_delta

        strain_dict ={}
        #Finally, grab each stress for every 0.01 value of strain
        strain_list = [round(num, 4).item() for num in np.linspace(0, 1, 101)]
        for number in strain_list:
            strain_diff = strain_loading.copy() - number
            closest_strain_idx = strain_diff.abs().idxmin()
            try:
                closest_strain_value = strain_loading.iloc[closest_strain_idx]
                closest_stress_value = stress_loading.iloc[closest_strain_idx]
                strain_dict[number] = closest_stress_value
            except IndexError:
                #If there's an index error, just pass on and hope for the best
                # strain_dict[number] = closest_stress_value
                strain_dict[number] = 'nan'
        
        
        #Whatever the last 'strain_dict' to be processed was, pass that as the last replicate stress list
        replicate_dict['last_replicate_stress_dict'] = strain_dict

            
        #Create output DataFrame for this replicate
        output_df = pd.DataFrame({
            'strain_data_loading': strain_loading,
            'strain_data_unloading': strain_unloading,
            'stress_data_loading': stress_loading,
            'stress_data_unloading': stress_unloading,
            'extension_data_loading': extension_loading,
            'extension_data_unloading': extension_unloading,
            'strain_lockup': strain_lockup.item(),      #convert to python native
            'extension_lockup': ext_lockup.item(),      #convert to python native
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
    name_split_list_under = filename.split('_')
    name_split_list_hyphen = filename.split('-')

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
    namerow_dict = check_for_namerow(mech_data_filepath)
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
        extension_name = this_parsing_dict['extension_column_name']
        stress_name_exclusions = this_parsing_dict['stress_exclusion_terms']
        strain_name_exclusions = this_parsing_dict['strain_exclusion_terms']
        extension_name_exclusions = this_parsing_dict['extension_exclusion_terms']

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
                if stress_name.lower() in cell_text.lower():
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in stress_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        stress_col_idx = col_idx
                        stress_col_name = cell_text
                if strain_name.lower() in cell_text.lower():
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in strain_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        strain_col_idx = col_idx
                        strain_col_name = cell_text

                if extension_name.lower() in cell_text.lower():
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in extension_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        extension_col_idx = col_idx
                        extension_col_name = cell_text

            #Finally, try to open the data
            try:
                data_df = pd.read_csv(mech_data_filepath, usecols=[stress_col_idx, strain_col_idx, extension_col_idx], skiprows=(start_row_idx) )
            except UnboundLocalError:
                #If you don't find a 'stress_column_name' or 'strain_column_name', assume this isn't a good mechanical file
                file_output_dict['dataframe'] = pd.DataFrame({})
                file_output_dict['parse_status']['file_pandas_readable'] = False
                file_output_dict['parse_status']['namerow_good'] = False
                file_output_dict['parse_status']['namerow_idx'] = 0
                
                return file_output_dict

            #Make sure the strain is right-side up (i.e. positive values only), ignore lead-in if not at 0, and set minimum at 0
            try:
                strain_series = data_df[strain_col_name].copy()
            except KeyError:
                #if keyerror, parsing has failed and take a look at the data_df to see what's going on
                print("Strain column parse fail; check 'data_df' below:")
                print(f"\t {data_df.head(10)}")
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
            strain_series = strain_series - new_min  #reset series minimum to a more useful 0
              # re-assign strain series
            data_df[strain_col_name] = strain_series

            #Make sure the extension is right-side up (i.e. positive values only), ignore lead-in if not at 0, and set minimum at 0
            try:
                extension_series = data_df[extension_col_name].copy()
            except KeyError:
                #if keyerror, parsing has failed and take a look at the data_df to see what's going on
                print("Extension column parse fail; check 'data_df' below:")
                print(f"\t {data_df.head(10)}")
              # get info about strain
            pos_ext_sum = len(extension_series[extension_series>0])
            neg_ext_sum = len(extension_series[extension_series<0])
            this_min_ext = extension_series.min()
              # flip strain if negative
            if neg_ext_sum > pos_ext_sum:
                # adj_ext_series = extension_series -this_min_ext
                # adj_ext_series = abs(this_min_ext)-adj_ext_series
                adj_ext_series = extension_series * -1
                extension_series = adj_ext_series  #reset series
              # check if a lead-in is artificially skewing 0
            data_start_idx = round(len(extension_series)*0.25)  #ignore first part
            new_min = extension_series[data_start_idx::].min()  #find real '0' without lead-in garbage
            strain_series = strain_series - new_min  #reset series minimum to a more useful 0
              # re-assign strain series
            data_df[extension_col_name] = extension_series

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
    namerow_dict = check_for_namerow(mech_data_filepath)
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
        extension_name = this_parsing_dict['exclusion_column_name']
        stress_name_exclusions = this_parsing_dict['stress_exclusion_terms']
        strain_name_exclusions = this_parsing_dict['strain_exclusion_terms']
        extension_name_exclusions = this_parsing_dict['extension_exclusion_terms']

        #Handle XLSX
        if mech_data_filepath.lower().endswith('.xlsx'):
            #Open the CSV and find the column name row
            excel_df = pd.read_excel(mech_data_filepath, sheet_name= good_sheetname)
            rows = excel_df.values.tolist()
              # assign header row to column names; if these are superseded by cell text, 'header_row' will be updated
              # i.e. if 'rows' DataFrame has garbage column names because the first row doesn't contain actual column names, the for loop below will find a better guess
            header_row = list(excel_df.columns)
            column_names_good = True
            for row_num, row in enumerate(rows):
                #Decrement 'start_row_idx' to account to switch in index start values
                if row_num == (start_row_idx-1):
                    header_row = row
                    column_names_good = False
            
            #Walk through header row and find which column index has the stress and strain columns
            for col_idx, cell_text in enumerate(header_row):
                if stress_name.lower() in str(cell_text).lower():
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in stress_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        stress_col_idx = col_idx
                        stress_col_name = cell_text
                if strain_name.lower() in str(cell_text).lower():
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in strain_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        strain_col_idx = col_idx
                        strain_col_name = cell_text
                if extension_name.lower() in cell_text.lower():
                    #Make sure no exlcusion terms are in the cell text
                    exclusion_sum = sum([1 for term in extension_name_exclusions if term.lower() in cell_text.lower()])
                    if exclusion_sum == 0:
                        extension_col_idx = col_idx
                        extension_col_name = cell_text

            #Finally, open the data
            if column_names_good:
                data_df = excel_df[[stress_col_name, strain_col_name]]
              # if the column names were found somewhere other than the first row, try and open only the stress/strain columns starting at that row
              # NOTE: this parsing can be janky on Excel files with specific formatting (i.e. automated instrument output reports)
            else:
                data_df = pd.read_excel(mech_data_filepath, usecols=[stress_col_idx, strain_col_idx], skiprows=(start_row_idx), sheet_name= good_sheetname)

            #Make sure the strain is right-side up (i.e. positive values only), ignore lead-in if not at 0, and set minimum at 0
            try:
                strain_series = data_df[strain_col_name].copy()
            except KeyError:
                #if keyerror, parsing has failed and take a look at the data_df to see what's going on
                print("Strain column parse fail; check 'data_df' below:")
                print(f"\t {data_df.head(10)}")
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
            strain_series = strain_series - new_min  #reset series minimum to a more useful 0
              # re-assign strain series
            data_df[strain_col_name] = strain_series

            #Make sure the extension is right-side up (i.e. positive values only), ignore lead-in if not at 0, and set minimum at 0
            try:
                extension_series = data_df[extension_col_name].copy()
            except KeyError:
                #if keyerror, parsing has failed and take a look at the data_df to see what's going on
                print("Extension column parse fail; check 'data_df' below:")
                print(f"\t {data_df.head(10)}")
              # get info about strain
            pos_ext_sum = len(extension_series[extension_series>0])
            neg_ext_sum = len(extension_series[extension_series<0])
            this_min_ext = extension_series.min()
              # flip strain if negative
            if neg_ext_sum > pos_ext_sum:
                # adj_ext_series = extension_series -this_min_ext
                # adj_ext_series = abs(this_min_ext)-adj_ext_series
                adj_ext_series = extension_series * -1
                extension_series = adj_ext_series  #reset series
              # check if a lead-in is artificially skewing 0
            data_start_idx = round(len(extension_series)*0.25)  #ignore first part
            new_min = extension_series[data_start_idx::].min()  #find real '0' without lead-in garbage
            strain_series = strain_series - new_min  #reset series minimum to a more useful 0
              # re-assign strain series
            data_df[extension_col_name] = extension_series

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


def split_filename_for_printname_guessing(filename, logbook_df,
                                          printname_examples = None):
    """
    Description:
        Take a mechanical data filename (CSV or XLSX), clean it up for comparison, and compare against the logbook to find a printname.

    INPUT:
        'filename'              str filepath
        'printname_examples'    list, pd.Series, np.array; take an iterable and make an np.array of lowrcase strings out of it
    ACTION:
        -lorem
    OUTPUT:
        'printname_dict'    lorem
        
    """

    #Initialize variables
      # format examples as an array from whatever iterable is passed (should be short iterable so speed doesn't matter)
    if printname_examples:
        printname_examples = np.array([str(entry).lower() for entry in printname_examples])
        printname_sets = [set(str(entry).lower()) for entry in printname_examples]
    else:
        printname_examples = logbook_df['Name']
        printname_examples = np.array([str(entry).lower() for entry in printname_examples])
        printname_sets = [set(str(entry).lower()) for entry in printname_examples]
      # make sure 'filename' isn't a filepath and format for comparison
    filename = os.path.basename(filename)
    filename = filename.lower().replace(".csv","").replace(".xlsx","")
      # initalize the return dict
    printname_dict = {
        'full_input_name': filename,
        'best_guess_printname': None,
        'printname_parse_bool': False,
        'iteration_marker': '',
        'all_matches': [],
        'all_match_scores': [],
        'logbook_entry': None
        }
      # regular expression pattern definitions for better printname parsing
    datetime_pattern = r"^\d{8}$"
    operator_pattern = r"[a-zA-Z]{2,3}"
    datetime_string = ''
    operator_strings = []
    iterator_strings = []


    #Split the filename to match patterns
    hyphen_list = filename.split('-')
    underscore_list = filename.split('_')

    #Check for formatting to match underscore format in logbook
    if len(hyphen_list) > len(underscore_list):
        #Run through each part of the filename and try to find components of a printname
        for idx, name_part in enumerate(hyphen_list):
            if re.match(datetime_pattern, name_part):
                datetime_string = name_part
            elif re.match(operator_pattern, name_part):
                operator_strings.append(name_part)
            else:
                iterator_strings.append(name_part)
    #Explicitly check for more underscores in case both split lists are short (i.e. 0 or 1)
    elif len(hyphen_list) < len(underscore_list):
        #Run through each part of the filename and try to find components of a printname
        for idx, name_part in enumerate(underscore_list):
            if re.match(datetime_pattern, name_part):
                datetime_string = name_part
            elif re.match(operator_pattern, name_part):
                operator_strings.append(name_part)
            else:
                iterator_strings.append(name_part)

    #Put together combinations of printname and check logbook print names for a match
    # set is faster than alternatives, so do a quick first pass
    best_name_guesses = []
    max_score = 0
    potential_matches = []
    match_scores = []
    for operator_string in operator_strings:
        for iterator_string in iterator_strings:
            this_name_guess = f"{datetime_string}_{operator_string}_{iterator_string}"
            this_name_set = set(this_name_guess)
            for idx, print_set in enumerate(printname_sets):
                bool_check = bool(
                    (len(this_name_set.intersection(print_set))>5) and
                    (datetime_string in printname_examples[idx])
                    )
                if bool_check:
                    this_score = len(this_name_set.intersection(print_set))/len(this_name_set)
                    potential_matches.append(printname_examples[idx])
                    match_scores.append(this_score)
                    if this_score >= max_score:
                        best_name_guesses.append(this_name_guess)
                        max_score = this_score

    #Run through matches and try to find the best one (second pass)
    # use entropy here because the list is much shorter
    # NOTE: this is portional entropy (i.e. it's non-negative and not zero)
    min_diff_score = 1   #hard-coded minimum score based on testing; anything higher than this is usually garbage
    best_name_guess = ''
    ent_matches = []
    ent_scores = []
    for print_name_guess in best_name_guesses:
        #Get portional entropies (values are 0-pos_inf theoretically, i.e. non-negative) of the target string ('print_name_guess')
        #   relative to the portional entropy of the matches.
        #   Use the match string characters to calculate entropy, otherwise two close matches could have the same entropy.
        for match in potential_matches:
            #Get character sets
            guess_set = set(print_name_guess)
            match_set = set(match)
            #Get letter counts and entropy for match string
            match_set_cnts = [match.count(character) for character in guess_set]
            match_set_probs = [cnt/sum(match_set_cnts) for cnt in match_set_cnts]
              # calculate the log if cnt isn't 0, return 0 if the cnt is <1 and math.log has a domain error
            match_set_logs = [math.log(cnt) if (cnt>0) else 0 for cnt in match_set_cnts ]
            match_port_ent = sum([prob*log for prob, log in zip(match_set_probs, match_set_logs)])
            #Get letter counts from the match letter set for the target string
            guess_set_cnts = [print_name_guess.count(character) for character in guess_set]
            guess_set_probs = [cnt/sum(guess_set_cnts) for cnt in guess_set_cnts]
            guess_set_logs = [math.log(cnt) if (cnt>0) else 0 for cnt in guess_set_cnts]
            guess_port_ent = sum([prob*log for prob, log in zip(guess_set_probs, guess_set_logs)])
            #Calculate a separate score metric based on how many of the target letter cnts are the same in the 'match' string
            cnt_diffs = [round((match_cnt/guess_cnt), 5) for match_cnt, guess_cnt in zip(match_set_cnts, guess_set_cnts)]
            #Compare entropies 
            relative_ent_diff = 1-(guess_port_ent/match_port_ent)
              # ideally 'cnt_diffs' is [1,1,1, ...] for a perfect match, so subtract length of 'guess_set' to get ideal match value of 0
              # deviation for each character from ideal could be >1 or <1, so 'relative_cnt_diff' magnitude is important
              # sign of 'cnt_diff_total' is a function of mostly string length; + if 'match' has a lot more characters but is a match, - if 'match' is shorter than expected
            relative_cnt_diff = sum(cnt_diffs)-len(guess_set)
            this_score = relative_cnt_diff
            #If score is better than threshold or previous best, add to everything important
            #   Adding in a buffer to catch close but off cases (i.e. the '1.2' value means close +/- 20%)
            if abs(this_score) <= (min_diff_score*2):
                #Save the scores and matches
                ent_matches.append(match)
                ent_scores.append(this_score)
                #Set best guess
                min_diff_score = abs(this_score)
                best_name_guess = print_name_guess

    #Set final return dict values as appropriate
    printname_dict['all_matches'] = ent_matches
    printname_dict['all_match_scores'] = ent_scores
    printname_dict['best_guess_printname'] = best_name_guess  #initial guess; return this guess if logbook lookup fails
    
    #Final printname guess and logbook lookup (third pass through printnames)
    if len(ent_matches) >0:
        #Take 'best_name_guess' printname, split, and find datacode index
        split_guess_name_list = best_name_guess.split('_')
        
        #Final format check before making assumptions about which data is in which spot
        datetime_idx = -1
        for idx, name_part in enumerate(split_guess_name_list):
            if re.match(datetime_pattern, name_part):
                datetime_string = name_part
                datetime_idx = idx
        clean_printname_iterator = ''
        iterator_junk_list = []
        if datetime_idx == 0:
            printname_iterator = split_guess_name_list[2]
            for idx, character in enumerate(printname_iterator):
                #look for last iterator format of "00"
                #   i.e. for '20260101_acp_01b', '_01' should be parsed as '01'
                if character.isdigit() and (idx <2):
                    clean_printname_iterator = clean_printname_iterator + str(character)
                else:
                    iterator_junk_list.append(character)
            #Finalize guess name
            clean_printname_guess = f"{split_guess_name_list[0]}_{split_guess_name_list[1]}_{clean_printname_iterator}"
            printname_dict['best_guess_printname'] = clean_printname_guess

            #Finalize logbook lookup
            try:
                if clean_printname_guess in printname_examples:
                    logbook_df['Name'] = logbook_df['Name'].str.lower()
                    printname_dict['logbook_entry'] = logbook_df[logbook_df['Name']==clean_printname_guess]
                printname_dict['printname_parse_bool'] = True
            except:
                pass

            #Filename remainder following logbook lookup
            remaining_junk = filename.replace(f"{split_guess_name_list[0]}", "")
            remaining_junk = remaining_junk.replace(f"{split_guess_name_list[1]}", "")
            remaining_junk = remaining_junk.replace(f"{clean_printname_iterator}", "")
            remaining_junk = remaining_junk.replace("-", "_")
            iterator_junk_list = [part for part in remaining_junk.split("_") if (len(part)>0)]

        #Finish up mechanical testing iterator guessing
        if len(iterator_junk_list) == 1:
            printname_dict['iteration_marker'] = iterator_junk_list[0]
          # if there's still junk left at this point, give up, assume the mechanical testing iteration marker is complicated, and just report it
        else:
            junk_string = ''
            for part in iterator_junk_list:
                junk_string = junk_string + '_' + str(character)
            printname_dict['iteration_marker'] = junk_string

    return printname_dict



