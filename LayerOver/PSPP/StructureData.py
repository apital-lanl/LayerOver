"""
Copyright 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2025-11-25
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Interact with and analyze data and predictions based on structure and materials (broadly defined). Interface library between 1) structure parameters and 2) property/prediction data.

"""
#Import libraries and modules
  # general Python libraries
import json
from tkinter import Tk, filedialog
from pathlib import Path
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import os
import shutil
import errno
  # LayerOver-specific modules

#######################################################################################################################
#####  Generic specification files  ###################################################################################
#######################################################################################################################

generic_materials_dict = {
        'll50': {
            'report_name':'',
            'density': 0,
            'viscosity': 0
            },
        'll60': {},
        'SE1700': {},
        'PDMS': {},
        'default': {
            }
        }

generic_ideality_dict = {
    'thickness': {},
    'density': {},
    'compression': {},
    }

#######################################################################################################################
#####  Generic Functions  #############################################################################################
#######################################################################################################################


def load_data_jsons(target_directory=None):
    '''
    Description: Get available catalogs of available data from configuration-defined filepaths to JSON files with summary info.
        If that lookup doesn't work because there is no configuration file, look for the JSON files directly.
        If no JSON summaries can be found, prompt user to select 1) a logbook file and 2) a directory to walk for mechanical data.

    INPUT:
      (optional)
        'target_directory'  str, filepath; if a specific directory is desired (for project control, for example), specify it.
                            Otherwise a standard root path and structure will be found or created.
    ACTIONS:
        - Look for LayerOver root
            - If none found, look for JSON config standards in User root
        - If JSON configs not found, create a LayerOver root
            - Create folders and template JSON configs
            - Prompt user to select desired logbook file
            - Prompot user to select a desired directory with logbook-associated mechanical data
        - If all else fails, raise error dialog for user that LayerOver root exists but is effectively empty
    OUTPUT:
        n/a
    '''
    
    pass


def find_data_by_structure_id(structure_name):

    #Attempt to find LayerOver root

        # if no root is found, create root


    #Load structure_dict from root

    pass


def parse_material_note(raw_material_note):
    """
    Description:
        lorem

    INPUT:
        'raw_material_note'         lorem
    ACTION:
        -lorem
    OUTPUT:
        'parsed_material_string'    materials from dict above
    """
    #Initialize variables
