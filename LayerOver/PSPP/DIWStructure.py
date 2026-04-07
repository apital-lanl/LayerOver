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

Description: Module to parse, homogenize, store, and compare DIW strucuture codes and name handling.

"""

#Import standard libraries
import numpy as np
import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from tkinter import Tk, filedialog
#Import LayerOver components
from LayerOver.PSPP.Materials import parse_material_note

#######################################################################################################################
#####  Template Variables  ############################################################################################
#######################################################################################################################

blank_diw_structure_dict = {
    'metadata': {
        'unique_structure_name': None,
        'print_name': None,
        'structure': None,
        'nozzle_size_um': None,
        'pitch_offset': None,
        'syringe-material': None,
        'project': None,
        'machine_name': None,
        'layerup_file': None,
        'version_number': None,
        'notes': None,
        'mech_data_note': None,
        'mech_data_filepaths': None,
        'keyence_data_note': None,
        'punch_diameter': None,
        'mass_g': None,
        'thickness_mm': None,
        'density_g/cc': None
        },
    'part_structure': None,
    'part_skin_nozzle_size': None,
    'part_layer_nozzle_size': None,
    'part_angular_offset': None,
    'part_lateral_offset': None,
    'part_material': '',
    'part_pitch': None,
    'number_of_layers': None,
    'layer_strand_extrusion': None,
    'layer_strand_diameter': None,
    'layer_types': None,
    'layer_type_modifiers': None,
    'layer_points': None,
    'layer_steps': None,
    'layer_angles': None,
    'layer_lateral_offsets': None,
    'layer_materials': None,
    'layer_pitches': None
    }

default_diw_structure_dict = {
    'metadata': {
        'unique_structure_name':'',
        'print_name': '',
        'structure': '',
        'nozzle_size_um': '',
        'pitch_offset': 0,
        'syringe-material': '',
        'project': '',
        'machine_name':'',
        'layerup_file': '',
        'version_number': '',
        'notes':'',
        'mech_data_note': '',
        'mech_data_filepaths': [],
        'keyence_data_note': '',
        'punch_diameter': '',
        'mass_g': 0,
        'thickness_mm': 0,
        'density_g/cc': 0
        },
    'part_structure':'',
    'part_skin_nozzle_size': 0,
    'part_layer_nozzle_size': 0,
    'part_angular_offset': 0,
    'part_lateral_offset': 0,
    'part_material': '',
    'part_pitch': 0,
    'number_of_layers': 0,
    'layer_strand_extrusion': 'constant',
    'layer_strand_diameter': [],
    'layer_types': [],
    'layer_type_modifiers': [],
    'layer_points': [],
    'layer_steps': [],
    'layer_angles': [],
    'layer_lateral_offsets': [],
    'layer_materials': [],
    'layer_pitches': []
    }

'''
metadata                    Metadata from logbook; if more than one entry exists for the given structure, lists are passed. Othewise str, int, or bool
part_structure              Generic 'S-code' structure; can be complete or short; parsed as string
part_skin_nozzle_size
part_layer_nozzle_size
part_angular_offset
part_lateral_offset
part_material
part_pitch
number_of_layers            number of layers in the part
layer_strand_extrusion      str ('constant', 'variable') or list (str for each layer ('constant', 'variable'), discreet values)
layer_strand_diameter       list; parsed based on 'layer_strand_extrusion' strand type
layer_types                 list of str for each layer ('helicoidal', 'spiral', 'maze')
    'helicoidal'             - generic layer for parallel strands; includes 'skin' layers; parsed as linear movements unless a 'layer_type_modifier' is used
    'spiral'                 - special layer type, parsed as an arc
    'maze'                   - generic flag for grid-based geometries
layer_type_modifiers        str ('None', 'Perturbed')
    'Perturbed'              - stochasticity added to the otherwise linear path
    'None'                   - just printed coordinate-to-coordinate with no perturbation or stochasticity added
layer_points                list of lists or numpy.array for each layer; initialized as empty list
layer_steps                 list of int/float; height of each layer from the substrate
layer_angles                list of int/float; angular offset for each layer relative to a global '0'; 
    eg. [0, 40, 80] would be a 40 degree offset for each subsequent layer to the previous
layer_lateral_offsets
layer_materials
layer_pitches
'''

#DIW material-to-opacity
diw_material_opacity_dict = {
    'default':{
        'opacity': 0.1,
        'per_length': 1,
        'length_unit': 'mm'},
    'll50': {
        'opacity': 0.1,
        'per_length': 1,
        'length_unit': 'mm'},
    'll60': {
        'opacity': 0.1,
        'per_length': 1,
        'length_unit': 'mm'},
    'SE1700': {
        'opacity': 0.1,
        'per_length': 1,
        'length_unit': 'mm'},
    'PDMS': {
        'opacity': 0.1,
        'per_length': 1,
        'length_unit': 'mm'},
    }

#str literals for structure
diw_structure_codes = {
    "s": {
        'name': 'skin',
        'parse_flag': 'helicoidal',
        },
    "h": {
        'name': 'helicoidal',
        'parse_flag': 'helicoidal',
        },
    "m": {
        'name': 'maze',
        'parse_flag': 'grid',
        }
    }

#Pulled straight from logbook entries
blank_diw_logbook_row_dict = {
    'print_name':'',
    'structure':'',
    'strand_diameter':'',
    'angle_of_rotation':'',
    'lateral_offset':'',
    'pitch':'',
    'syringe-material':'',
    'project':'',
    'thickness': '',
    'density': ''
    }

standard_logbook_columnnames = [
    'Name',	
    'Structure',
    r'Strand Diameter, nominal (skin/heli)',
    'Angle of Rotation (deg)',
    'Lateral Offset (um)',	            #encoding errors are a pain, so assume "um" always means micro-meters (10^-6 m)
    'Pitch (um)',	
    'Pitch Layer List',	
    r'Syringe/Material',	
    'Project',	
    'Machine Name',	
    'LayerUp File',	
    'Version #',
    'Notes',	
    'Mechanical Data? (initals)',	
    'Keyence? (initials)',	
    'Punch Diameter',	
    'Mass (g)',	
    'Thickness (Checkline) (mm)',	
    'Thickness (Confocal) (mm)',
    'Thickness (Fancy KCNSC) (mm)',	
    r'Density (g/cc)',	
    'Thickness (Additional) (mm)',	
    r'Thickness/ Density Initials',	
    'Humidity',	
    'Column1',
    'Strand Diameter, Skin',
    'Strand Diameter, Layer'
    ]

blank_unique_structure_dict = {
    'part_structure': None,
    'layer_strand_diameter': None,
    'layer_types': None,
    'layer_angles': None,
    'layer_lateral_offsets': None,
    'layer_materials': None,
    'layer_pitches': None
    }

#The following list items are additional fields that may be appended to logbook entries based on various operations
#  Possible operations: mechanical data import, opacity/volumetric prediction from structure, analysis of PSPP
additional_structure_fields = [
    'index',  #unique structure ID name from opacity/volumetric prediction
    'vox_1_average',
    'vox_1_sum',
    'vox_1_density',
    'vox_2_average',
    'vox_2_sum',
    'vox_2_density',
    'vox_3_average',
    'vox_3_sum',
    'vox_3_density',
    'full_voxel_dimensions',
    'array_side_length_mm',
    'vox_resolution_micron'
    ]


#######################################################################################################################
#####  Generic Functions  #############################################################################################
#######################################################################################################################




#######################################################################################################################
#####  DIW-specific Functions  ########################################################################################
#######################################################################################################################


def structure_dict_from_logbook_row(logbook_row):
    '''
    Description: Fill in a skeleton 'blank_diw_structure_dict' with values from a logbook DataFrame row.
        If no values are found, fill with 'default_diw_structure_dict' values. 'None' entries are failed entries. 
    '''

    #Initialize variables
    structure_dict = blank_diw_structure_dict.copy()
      #Make sure types are correct; some values simply re-set to allow for future re-typing as required
    logbook_row['Structure'] = logbook_row['Structure']
    logbook_row['Strand Diameter, Skin'] = logbook_row['Strand Diameter, Skin'].astype(float)
    logbook_row['Strand Diameter, Layer'] = logbook_row['Strand Diameter, Layer'].astype(float)
    logbook_row['Angle of Rotation (deg)'] = logbook_row['Angle of Rotation (deg)'].astype(float)
    logbook_row['Lateral Offset (µm)'] = logbook_row['Lateral Offset (µm)'].astype(float)
    logbook_row['Name'] = logbook_row['Name']
    logbook_row['Pitch (µm)'] = logbook_row['Pitch (µm)'].astype(float)
    logbook_row['Machine Name'] = logbook_row['Machine Name']
    logbook_row['LayerUp File'] = logbook_row['LayerUp File']
    logbook_row['Version #'] = logbook_row['LayerUp File']
    logbook_row['Notes'] = logbook_row['Notes']
    logbook_row['Mechanical Data? (initals)'] = logbook_row['Mechanical Data? (initals)']
    logbook_row['Punch Diameter'] = logbook_row['Punch Diameter']
    logbook_row['Mass (g)'] = logbook_row['Mass (g)'].astype(float)
    logbook_row['Thickness (Checkline) (mm)'] = logbook_row['Thickness (Checkline) (mm)'].astype(float)
    logbook_row['Thickness (Confocal) (mm)'] = logbook_row['Thickness (Confocal) (mm)'].astype(float)
    logbook_row['Thickness (Fancy KCNSC) (mm)'] = logbook_row['Thickness (Fancy KCNSC) (mm)'].astype(float)
    logbook_row['Density (g/cc)'] = logbook_row['Density (g/cc)'].astype(float)
    logbook_row['Thickness (Additional) (mm)'] = logbook_row['Thickness (Additional) (mm)'].astype(float)
    logbook_row['Thickness/ Density Initials'] = logbook_row['Thickness/ Density Initials'].astype(str)
    logbook_row['Humidity'] = logbook_row['Humidity']
    logbook_row['Column1'] = logbook_row['Column1']
    logbook_row['Pitch Layer List'] = logbook_row['Pitch Layer List'].astype(str)

    #Add direct metadata and other fields
    structure_dict['part_structure'] = logbook_row['Structure'].values[0]
    structure_dict['part_skin_nozzle_size'] = logbook_row['Strand Diameter, Skin'].values[0]
    structure_dict['part_layer_nozzle_size'] = logbook_row['Strand Diameter, Layer'].values[0]
    structure_dict['part_angular_offset'] = logbook_row['Angle of Rotation (deg)'].values[0]
    structure_dict['part_lateral_offset'] = logbook_row['Lateral Offset (µm)'].values[0]
    structure_dict['metadata']['print_name'] = logbook_row['Name'].values[0]
    structure_dict['metadata']['pitch_offset'] = logbook_row['Pitch (µm)'].values[0]
    structure_dict['metadata']['machine_name'] = logbook_row['Machine Name'].values[0]
    structure_dict['metadata']['layerup_file'] = logbook_row['LayerUp File'].values[0]
    structure_dict['metadata']['version_number'] = logbook_row['Version #'].values[0]
    structure_dict['metadata']['notes'] = logbook_row['Notes'].values[0]
    structure_dict['metadata']['mech_data_note'] = logbook_row['Mechanical Data? (initals)'].values[0]
    structure_dict['metadata']['punch_diameter'] = logbook_row['Punch Diameter'].values[0]
    structure_dict['metadata']['mass_g'] = logbook_row['Mass (g)'].values[0]
    structure_dict['metadata']['thickness_mm'] = logbook_row['Thickness (Checkline) (mm)'].values[0]
    # structure_dict['metadata'][''] = logbook_row['Thickness (Confocal) (mm)']
    # structure_dict['metadata'][''] = logbook_row['Thickness (Fancy KCNSC) (mm)']
    structure_dict['metadata']['thickness_mm'] = logbook_row['Density (g/cc)'].values[0]
    # structure_dict['metadata'][''] = logbook_row['Thickness (Additional) (mm)']
    # structure_dict['metadata'][''] = logbook_row['Thickness/ Density Initials']
    # structure_dict['metadata'][''] = logbook_row['Humidity']
    # structure_dict['metadata'][''] = logbook_row['Column1']
    structure_dict['layer_pitches'] = logbook_row['Pitch Layer List'].values[0]

    #Pull variables that require interpretation
      # parse structure code to get list of layer types and number of layers  
    layer_list = parse_structure(logbook_row['Structure'].values[0])
    number_of_layers = len(layer_list)
    structure_dict['layer_types'] = layer_list
    structure_dict['number_of_layers'] = number_of_layers
      # material
    raw_logbook_material_string = logbook_row['Syringe/Material'].values[0]
    parsed_material_note = parse_material_note(raw_logbook_material_string)
    split_material_note = parsed_material_note.split('_')
    base_material = split_material_note[0]
    layer_material_list = []
    if 'homogenous' in parsed_material_note.lower():
        structure_dict['metadata']['part_material'] = parsed_material_note
        for i in range(number_of_layers):
            layer_material_list.append(base_material)
    else:
        #TODO: add use case for parsing if multiple materials are used in the print
        structure_dict['metadata']['part_material'] = parsed_material_note
        for i in range(number_of_layers):
            pass
    structure_dict['layer_materials'] = layer_material_list   
      # parse the pitch 
      #TODO: add support for parsing the 'raw_pitch_string' 
    raw_pitch_string = structure_dict['metadata']['pitch_offset']
    pitch_entry = logbook_row['Pitch Layer List'].values[0]
    if type(pitch_entry) == str:
        pitch_list = []
        list_split = pitch_entry.split(',')
        if len(list_split) >1:
            pitch_list = [float(part.strip()) for part in list_split]
        else:
            for i, layer_type in zip(range(number_of_layers), structure_dict['layer_types']):
                if layer_type == 'skin':
                    pitch_list.append(structure_dict['part_skin_nozzle_size'])
                else:
                    pitch_list.append(float(pitch_entry))
    elif type(pitch_entry) == float:
        pitch_list = []
        for i, layer_type in zip(range(number_of_layers), structure_dict['layer_types']):
            if layer_type == 'skin':
                pitch_list.append(structure_dict['part_skin_nozzle_size'])
            else:
                pitch_list.append(float(pitch_entry))
    else:
        pitch_list = pitch_entry
    structure_dict['layer_pitches'] = pitch_list
      # layer strand diameters
    skin_diam = structure_dict['part_skin_nozzle_size']
    layer_diam = structure_dict['part_layer_nozzle_size']
    strand_diam_list = []
    for layer_type in layer_list:
        if layer_type.lower() == 'skin':
            strand_diam_list.append(skin_diam)
        else:
            strand_diam_list.append(layer_diam)
    structure_dict['layer_strand_diameter'] = strand_diam_list
      # lateral offsets
      #TODO: add support for variable offsets by layer from logbook files
    offset = structure_dict['part_lateral_offset']
    offset_list = []
    if type(offset)== np.float64:
        for i in range(number_of_layers):
            offset_list.append(offset)
    elif type(offset) == float:
        for i in range(number_of_layers):
            offset_list.append(offset)
    elif type(offset) == str:
        try:
            for i in range(number_of_layers):
                offset_list.append(offset)
        except:
            print(f"str fail; Defaulting to 0 lateral offset for entry {offset}")
            for i in range(number_of_layers):
                offset_list.append(0)
    else:
        print(f"type fail; Defaulting to 0 lateral offset for entry {offset}: type {type(offset)}")
        for i in range(number_of_layers):
            offset_list.append(0)
    structure_dict['layer_lateral_offsets'] = offset_list       
      # angular offsets
    angle = structure_dict['part_angular_offset']
    angle_list = []
    if type(angle)== float:
        for i in range(number_of_layers):
            if i == 0:
                angle_list.append(0)
            else:
                this_angle = (angle * i) % 360
                angle_list.append(this_angle)
    elif type(angle) == str:
        try:
            if angle.isnumeric():
                angle = float(angle)
                for i in range(number_of_layers):
                    if i == 0:
                        angle_list.append(0)
                    else:
                        this_angle = (angle * i) % 360
                        angle_list.append(this_angle)
        except:
            print(f"Failure in angular offset for angle {angle} of type {type(angle)}")
    elif (type(angle)== np.float64):
        for i in range(number_of_layers):
            if i == 0:
                angle_list.append(0)
            else:
                this_angle = (angle * i) % 360
                angle_list.append(this_angle)
    else:
        print(f"Failure in angular offset for angle {angle} of type {type(angle)}")
    structure_dict['layer_angles'] = angle_list

    #TODO: add better functionality to these entries
    structure_dict['layer_strand_extrusion'] = [None for i in range(number_of_layers)]

    return structure_dict


def structure_dict_from_param_kwargs(structure_code = None,
                                     nozzle_size = None,
                                     skin_nozzle_size = None,
                                     layer_nozzle_size = None,
                                     ):
    '''
    Description: Given a 'blank_diw_param_dict'-like dictionary, generate a skeleton 'blank_structure_dict'.
    '''
    #Attempt to pull 

    pass


def parse_structure(structure_string, flag = ''):
    '''
    Description: Take a short, long, or complete structure string and generate a list of it's layer codes.

    INPUT:
        'structure_string'  str or stringlike object
    ACTION:
        -Split string into layer components
    OUTPUT:
        'layer_list'        list; each item in list is a layer type ('skin', 'hellicoidal, 'maze')

    TODO:
        -Add parsing support for other structure strings. Currently only supports "S8H" type strings
    '''
    
    #Condition variables
    structure_string = str(structure_string)

    #Flag option for changing parsing behavior 
    if flag == '':
        underscore_check = bool( len(structure_string.split("_")) > 1)
        curly_bracket_check = bool( len(structure_string.split("{")) > 1)
        flag = 'default'

    #Split the structure string character-by-character
    initial_length = len(structure_string)
    structure_parts = []
    latest_part = ''
    if flag == 'default':
        for idx, character in enumerate(structure_string):
            #add any numeric character to the string
            if character.isnumeric():
                latest_part = latest_part + character
            elif character != '':
                if latest_part == '':
                    structure_parts.append(character)
                else:
                    latest_part = latest_part + character
                    structure_parts.append(latest_part)
                    latest_part = ''

    #Now that the structure is split out, generate a list with the layer types
    layer_list = []
    for part in structure_parts:
        #If part is only one character, let's hope it's one of the options below
        if len(part) ==1:
            if part.lower() == 's':
                layer_list.append('skin')
            if part.lower() == 'h':
                layer_list.append('helicoidal')
            if part.lower() == 'm':
                layer_list.append('maze')
        #Otherwise, if numerics are involved then parse to get number of layers and add that number to the return array
        else:
            layer_number_string = ''
            type_string = ''
            for character in part:
                if character.isnumeric():
                    layer_number_string = layer_number_string + character
                else:
                    if type_string == '':
                        type_string = character.lower()
                    else:
                        print(f"Failure to parse string '{part}': too many layer type designators (ie. 's', 'h', 'm'")
            number_of_layers = layer_number_string
            if type(number_of_layers) == str:
                if number_of_layers == "":
                    number_of_layers = 1
                else:
                    number_of_layers = int(number_of_layers)
            for i in range(number_of_layers):
                if type_string == 's':
                    layer_list.append('skin')
                if type_string == 'h':
                    layer_list.append('helicoidal')
                if type_string == 'm':
                    layer_list.append('maze')

    return layer_list


def match_structure_to_unique_name(structure_dict, unique_structure_dict):
    '''
    Description: Take a set of structural parameters and return existing prints that are matches with an associated match score.

    OUTPUT:
        is_unique_bool          bool; was a match found or not
        unique_structure_id     str; unique structure ID for this run; not global
        unique_structure_dict   dict; all unique structures with keys of unique structure id (e.g. SHS-1, S5HS-9, etc.)

    '''
    #Initialize variables
    unique_names = list(unique_structure_dict.keys())
    this_structure = structure_dict['part_structure']
    is_unique_bool = False   #"Is the 'structure_dict' already in the 'unique_structure_dict'
      # initialize a list of structure fields to compare
    structure_fields = list(blank_unique_structure_dict.keys())
    
    #Create a dict with each structure's highest id number up to this point
    unique_structure_numbers = {}
    last_id_number = 0
    for name in unique_names:
        structure, structure_id = name.split('_')
        try:
            this_id_number = unique_structure_numbers[structure]
            if structure_id > this_id_number:
                unique_structure_numbers[structure] = structure_id
        except KeyError:
            unique_structure_numbers.update({structure: structure_id})

    
    #Compare the passed structure with all other structures
    match_found = False
    matching_structure = ''
    for known_structure_id in unique_names:
        comparison_structure_dict = unique_structure_dict[known_structure_id]
        matching_structure_bool = True
        #Compare each structure parameter set with all previous, unique structures
        for structure_field in structure_fields:
            comp_bool = comparison_structure_dict[structure_field] == structure_dict[structure_field]
            #boolean AND between 'matching' flag and each field's match boolean above
            matching_structure_bool &= comp_bool 
        #Match found
        if matching_structure_bool:
            match_found = True
            matching_structure = known_structure_id
            #return 'False', <match_structure>, <pass-through original dict>
            return is_unique_bool, known_structure_id, unique_structure_dict
    
    #No match found; update 'unique_structure_dict, find next index for unique structure id, flip 'is_unique_bool'
    if not match_found:
        #Create a blank structure dict to add to passed dict
        new_entry = blank_unique_structure_dict.copy()
        #Fill in each structure parameter field to be written to the passed dict below
        for structure_field in structure_fields:
            new_entry[structure_field] = structure_dict[structure_field]
        #Look for latest index number to increment and append for the passed dict key
        try:
            last_index = unique_structure_numbers[this_structure]
        #If 'KeyError', add as first iteration of that global structure
        except KeyError:
            last_index = 0
        #Create the new unique structure ID and add to passed structure dict
        new_index = int(last_index) +1
        unique_structure_numbers[this_structure]= new_index
        new_structure_id = this_structure+'_'+str(new_index)
        unique_structure_dict.update({new_structure_id:new_entry})
        is_unique_bool = True
        #return 'True', <new unique structure id>, <updated structure dict>
        return is_unique_bool, new_structure_id, unique_structure_dict

    #Catch-all case; shouldn't ever evaluate unless something goes extremely wrong
    return is_unique_bool, 'failure', unique_structure_dict


def structure_name_from_structure_dict(structure_dict):
    '''
    Directory:
        Lorem.
    '''
    
    #Initialize variables
    name_dict = {
        'generic_structure':'',
        'unique_structure_name': '',
        'default_structure_name': '', 
        }

    return name_dict


def open_wafflebook(waffledata_filepath,
                    preview_loaded_df = False):
    """
    Description:
        Open a '*WaffleData.csv' file output from a set of VolumetricPrediction runs.
    """

    #Initialize variables
    if os.path.exists(waffledata_filepath):
        pass
    else:
        print()
        print('Invalid filepath given for WaffleData; prompting for a good filepath.')
        root = Tk()
        waffledata_filepath = filedialog.askopenfilename(title = "Select a logbook 'WaffleData' file", \
                                                      filetypes = [('Waffle Structure Summary Files', '*WaffleData.csv')])
        root.destroy()

    waffledata_df = pd.read_csv(waffledata_filepath, encoding_errors= 'ignore')

    if preview_loaded_df:
        print(waffledata_df.head(10))

    return waffledata_df