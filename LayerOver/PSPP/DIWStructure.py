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
        'mech_data_flag': None,
        'keyence_data_flag': None,
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
        'mech_data_flag': False,
        'keyence_data_flag': False,
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

    #Add direct metadata and other fields
    structure_dict['part_structure'] = logbook_row['Structure']
    structure_dict['metadata'][''] = logbook_row['Strand Diameter, Skin']
    structure_dict['metadata'][''] = logbook_row['Strand Diameter, Layer']
    structure_dict['metadata'][''] = logbook_row['Angle of Rotation (deg)']
    structure_dict['metadata'][''] = logbook_row['Lateral Offset (µm)']
    structure_dict['metadata'][''] = logbook_row['Pitch (µm)']
    structure_dict['metadata'][''] = logbook_row['Syringe/Material']
    structure_dict['metadata'][''] = logbook_row['Project ']
    structure_dict['metadata'][''] = logbook_row['Machine Name']
    structure_dict['metadata'][''] = logbook_row['LayerUp File']
    structure_dict['metadata'][''] = logbook_row['Version #']
    structure_dict['metadata'][''] = logbook_row['Notes']
    structure_dict['metadata'][''] = logbook_row['Mechanical Data? (initals)']
    structure_dict['metadata'][''] = logbook_row['Punch Diameter']
    structure_dict['metadata'][''] = logbook_row['Mass (g)']
    structure_dict['metadata'][''] = logbook_row['Thickness (Checkline) (mm)']
    structure_dict['metadata'][''] = logbook_row['Thickness (Confocal) (mm)']
    structure_dict['metadata'][''] = logbook_row['Thickness (Fancy KCNSC) (mm)']
    structure_dict['metadata'][''] = logbook_row['Density (g/cc)']
    structure_dict['metadata'][''] = logbook_row['Thickness (Additional) (mm)']
    structure_dict['metadata'][''] = logbook_row['Thickness/ Density Initials']
    structure_dict['metadata'][''] = logbook_row['Humidity']
    structure_dict['metadata'][''] = logbook_row['Column1']
    structure_dict['layer_pitches'] = logbook_row['Pitch Layer List']

    #Pull variables
    

    #Condition row data for parsing and define terms  
    layer_list = parse_structure(row['Structure'])
    layer_length = len(layer_list)


    # default_diw_structure_dict = {
    # 'metadata': {
    #     'unique_structure_name':'',
    #     'print_name': '',
    #     'structure': '',
    #     'nozzle_size_um': '',
    #     'pitch_offset': 0,
    #     'syringe-material': '',
    #     'project': '',
    #     'machine_name':'',
    #     'layerup_file': '',
    #     'version_number': '',
    #     'notes':'',
    #     'mech_data_flag': False,
    #     'keyence_data_flag': False,
    #     'punch_diameter': '',
    #     'mass_g': 0,
    #     'thickness_mm': 0,
    #     'density_g/cc': 0
    #     },
    # 'part_structure':'',
    # 'part_skin_nozzle_size': 0,
    # 'part_layer_nozzle_size': 0,
    # 'part_angular_offset': 0,
    # 'part_lateral_offset': 0,
    # 'part_material': '',
    # 'part_pitch': 0,
    # 'number_of_layers': 0,
    # 'layer_strand_extrusion': 'constant',
    # 'layer_strand_diameter': [],
    # 'layer_types': [],
    # 'layer_type_modifiers': [],
    # 'layer_points': [],
    # 'layer_steps': [],
    # 'layer_angles': [],
    # 'layer_lateral_offsets': [],
    # 'layer_materials': [],
    # 'layer_pitches': []
    # }

    return 


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
                structure_parts.append(latest_part)
                latest_part = latest_part + character
            else:
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
            number_of_layers = int(layer_number_string)
            for i in range(number_of_layers):
                if type_string == 's':
                    layer_list.append('skin')
                if type_string == 'h':
                    layer_list.append('helicoidal')
                if type_string == 'm':
                    layer_list.append('maze')

    return layer_list


def match_structure_to_name(structure_dict):
    '''
    Description: Take a set of structural parameters and return existing prints that are matches with an associated match score.
    '''
    
    pass


def parse_logbook_row_for_structure(row_array):
    '''
    Directory:
        Lorem.
    '''
    
    #Initialize variables

    return structure_dict


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