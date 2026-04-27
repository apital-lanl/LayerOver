"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-04-23
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Higher-level functions for loading and implementing mechanical data analysis done with Analysis.MechanicalData.

"""

#Import libraries
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from tkinter import filedialog, Tk


  # import other LayerOver modules or module parts
from LayerOver.PSPP.DIWStructure import blank_diw_logbook_row_dict
from LayerOver.PSPP.DIWStructure import standard_logbook_columnnames
from LayerOver.PSPP.DIWStructure import parse_structure
from LayerOver.PSPP.DIWStructure import structure_dict_from_logbook_row


def open_mechanical_summary(mechsummary_filepath):
    """
    Description:
        Open a mechanical data summary file produced from Examples.Mechanical-Files_Directory_Process.py 

    INPUT:
        'mechsummary_filepath'  str (filepath); 
    ACTION:
        -lorem
    OUTPUT:
        'mechdata_df'        pd.DataFrame; 'cleaned' summary with standardized column data and column names
    """

    #Initialize variables
    mechdata_df = pd.DataFrame([])

    #Open the logbook
    if mechsummary_filepath.lower().endswith('.csv'):
        mechdata_df = pd.read_csv(mechsummary_filepath, encoding_errors= 'ignore')
        header_row = list(mechdata_df.columns)
    else:
        print(f"Check that filepath leads to a proper 'Summary of Directory Mechanical Data.csv' file")
        return mechdata_df

    #Add columns and format
    # mechdata_df['new_column'] = None
    # mechdata_df['new_column'] = mechdata_df['new_column'].astype('object')
    
    #Run through and directly add 
      # get the column names for columns to be cleaned
    for column_name in header_row:
        pass

    
    #Keeping the following in for reference, but no need to process these columns as much as for logbook
    # #Clean columns that tend to have garbage
    # logbook_df['Thickness (Checkline) (mm)'] = logbook_df['Thickness (Checkline) (mm)'].apply(lambda x: float(x) if (isinstance(x, (int, float))) else np.nan)

    # #Split out skin vs. layer nozzle size if different
    # logbook_df['Strand Diameter, Skin'] = [entry[0] if (len(entry)>1) else entry[0] for entry in logbook_df[diameter_column_name].apply(lambda s: str(s).split(r'/'))]
    # logbook_df['Strand Diameter, Skin'] = logbook_df['Strand Diameter, Skin'].astype(float)

    # #Clean the pitch column
    # #'Pitch Layer List' and 'layer_pitch_list'
    #   # if pitch column as an 'x', replace cell value with layer nozzle size times 'x' ammount (i.e. '1.25x' becomes "1.25 * layer_nozzle_size")
    # logbook_df[pitch_column_name] = logbook_df[pitch_column_name].astype(str)
    # x_mask = logbook_df[pitch_column_name].str.contains('x', case=False, na=False)
    # logbook_df.loc[x_mask, pitch_column_name] = (
    #     logbook_df.loc[x_mask, pitch_column_name].str.replace('x', '', case=False).astype(float) * 
    #     logbook_df.loc[x_mask, 'Strand Diameter, Layer']
    #     )
    # logbook_df.loc[x_mask, pitch_list_column_name] = logbook_df.loc[x_mask, pitch_column_name]

    return mechdata_df


def add_mechdata_to_logbook_df(logbook_df, mechdata_df):
    """
    Description:
        Add mechanical data from a mechdata_df (from open_mechanical_summary) to a logbook_df (from open_logbook) by matching on 'Structure Name' and 'Layer Name' columns. 
        This is done by creating new columns in the logbook_df for each column in the mechdata_df (except for 'Structure Name' and 'Layer Name') and filling those columns with the corresponding data from the mechdata_df where there are matches on 'Structure Name' and 'Layer Name'. 
    INPUT:
        'logbook_df'       pd.DataFrame; logbook dataframe to be updated with mechanical data
        'mechdata_df'      pd.DataFrame; mechanical data summary dataframe to pull data from)
    ACTION:
        -lorem
    OUTPUT:
        'logbook_df'       pd.DataFrame; updated logbook dataframe with new columns for mechanical data
    """

    #"Summary of Directory Mechanical Data.csv" columns as of 2026-04-23
    #   filename	filepath	filetype	immediate_parent_directory	parse_type	clean_filename	
    #   successful_columname_parse	file_namerow_index	pandas_readable	strain_assymptote	
    #   extension_assymptote	stress_at_0.01_strain	stress_at_0.02_strain	stress_at_0.03_strain	
    #   stress_at_0.04_strain	stress_at_0.05_strain	stress_at_0.06_strain	stress_at_0.07_strain	
    #   stress_at_0.08_strain	stress_at_0.09_strain	stress_at_0.1_strain	stress_at_0.2_strain	
    #   stress_at_0.3_strain	stress_at_0.4_strain	stress_at_0.5_strain	stress_at_0.6_strain	
    #   stress_at_0.7_strain	stress_at_0.8_strain	stress_at_0.85_strain	stress_at_0.9_strain	
    #   logbook_entry_found	mechanical_data_iteration	Name	Structure	Strand Diameter, nominal (Âµm) (skin/heli)	
    #   Angle of Rotation (deg)	Lateral Offset (Âµm)	Pitch (Âµm)	Layer 1 Height Multiplier	Layer 2+ Height Multiplier	
    #   layer_diameter_list	layer_angle_list	layer_lateral_list	layer_pitch_list	layer_height_modifier_list	
    #   Syringe/Material	Project	Machine Name	LayerUp File	Version #	
    #   Notes	Mechanical Data? (initals)	Keyence? (initials)	Punch Diameter	Mass (g)	
    #   Thickness (Checkline) (mm)	Thickness (Confocal) (mm)	Thickness (Fancy KCNSC) (mm)	
    #   Density (g/cc)	Thickness (Additional) (mm)	Thickness/ Density Initials	Humidity	
    #   Column1	Strand Diameter, Skin	Strand Diameter, Layer

    #Pull 'Name' column from mech data summary (direct logbook printnames)
    mechdata_printnames = mechdata_df['Name'].unique().tolist()
    logbook_df['mech_data_directory'] = None
    logbook_df['lower_name'] = logbook_df['Name'].apply(lambda s: str(s).lower())

    #Run through each unique printname and add data to logbook_df
    for idx, printname in enumerate(mechdata_printnames):
        #Pull mech data for this printname
        mechdata_printname_df = mechdata_df[mechdata_df['Name'] == printname]
        #Pull logbook rows with matching printname
        logbook_printname_mask = logbook_df['lower_name'] == printname
        try:
            #Pull directory mech data for printname was found in and add to logbook_df
            mechdata_filepath = mechdata_printname_df['filepath'].values[0]
            mechdata_dir = os.path.dirname(mechdata_filepath)
            logbook_df.loc[logbook_printname_mask, 'mech_data_directory'] = mechdata_dir
        except:
            logbook_df.loc[logbook_printname_mask, 'mech_data_directory'] = 'fail'
            print("Error in mech data filepath parsing (PSPP.DIWMechanicalProperties)")
            print()

        try:
            #Pull mechanical data replicates and average stress and extension data
            mean_strain_assymptote = mechdata_printname_df['strain_assymptote'].mean(numeric_only = True)
            mean_extension_assymptote = mechdata_printname_df['extension_assymptote'].mean(numeric_only = True)
            logbook_df.loc[logbook_printname_mask, 'mean_strain_assymptote'] = mean_strain_assymptote.astype(float)
            logbook_df.loc[logbook_printname_mask, 'mean_extension_assymptote'] = mean_extension_assymptote.astype(float)
            logbook_df.loc[logbook_printname_mask, 'n_mech_replicates'] = len(mechdata_printname_df)

        except:
            strain_means = mechdata_printname_df['strain_assymptote'].values
            ext_means = mechdata_printname_df['extension_assymptote'].values
            logbook_df[logbook_printname_mask]['mean_strain_assymptote'] = 'fail'
            logbook_df[logbook_printname_mask]['mean_extension_assymptote'] = 'fail'
            logbook_df[logbook_printname_mask]['n_mech_replicates'] = 'fail'
            print("Error in mech data replicate parsing (PSPP.DIWMechanicalProperties)")

    return logbook_df


