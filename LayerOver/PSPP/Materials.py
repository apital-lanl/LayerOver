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



#######################################################################################################################
#####  Generic specification data  ####################################################################################
#######################################################################################################################

generic_materials_dict = {
        'LL50': {
            
            },
        'LL60': {},
        'SE1700': {},
        'PDMS': {},
        'default': {
            'proper_name':'',
            'spec_type': 'broad',
            'density_mean': 0,
            'density_stdev': 0,
            'viscosity_mean': 0,
            'viscosity_stdev': 0,
            }
        }

    #units and standards
    #   'report_name'   str; name to use when key is found within string
    #   'density'       list of floats; 

#terms used in a material note and the 'material_name' they refer to
material_designations = {
    'default': {
        'material_name': 'defualt',
        '':'',},
    'll50': {
        'material_name': 'LL50',
        '':'',},
    'getl50': {
        'material_name': 'LL50, getter-filled',
        '':'',},
    'll60': {
        'material_name': 'LL60',
        '':'',},
    'tar': {
        'material_name': 'TAR',
        '':'',},
    'se1700': {
        'material_name': 'SE1700',
        '':'',},
    'lkr': {
        'material_name': 'LKR',
        '':'',},
    'pdms': {
        'material_name': 'PDMS',
        '':'',},
    'llama 50': {
        'material_name': 'LL50',
        '':'',},
    'llama 60': {
        'material_name': 'LL60',
        '':'',},
    }

#keys are terms found within a note about materials
material_modifiers = {
    'helical fibers': {
        'modifier_name': 'HelFib'
        },
    'default': {
        'modifier_name': 'default'
        },
    }



#######################################################################################################################
#####  Generic Functions  #############################################################################################
#######################################################################################################################


def parse_material_note(raw_material_note, 
                        layer_specification_notes = None):
    """
    Description:
        Parse logbook 'Syringe/Material' or related material note

    INPUT:
        'raw_material_note'         lorem
    ACTION:
        -lorem
    OUTPUT:
        'parsed_material_string'    materials from dict above; format "<basematerial>_<materialmodifiers_...>_<layer_specifics?>
    
    TODO:
        -Add parsing for layer-dependent material options
        
    """
    #Initialize variables
    raw_material_note = str(raw_material_note)
    parsed_material_string = ''

    #Get materials and modifiers
    material_list = []
    for key in list(material_designations.keys()):
        if key.lower() in raw_material_note.lower():
            material_name = material_designations[key]['material_name']
            if material_name not in material_list:
                material_list.append(material_name)
    modifier_list = []
    for key in list(material_modifiers.keys()):
        if key.lower() in raw_material_note.lower():
            material_mod = material_modifiers[key]['modifier_name']
            if material_mod not in modifier_list:
                modifier_list.append(material_mod)

    #Combine materials and modifiers
    base_material_string = ''
    for material in material_list:
        if base_material_string == '':
            base_material_string = material
        else:
            base_material_string = base_material_string + '-' + material
    base_modifier_string = ''
    for mod in modifier_list:
        base_modifier_string = base_modifier_string + '_' + mod

    #Add layer specifics as applicable
    if layer_specification_notes == None:
        parsed_material_string = base_material_string + base_modifier_string + '_homogenous'
    else:
        #TODO: add this use case to allow for different materials at different layers
        pass

    return parsed_material_string

