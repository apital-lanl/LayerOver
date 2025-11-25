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

Description: Module for handling conversions between ideal structure assumptions and physical measurements of as-printed structures.
    NOTE: Code below is for specific nomenclature and testing conventions. An attempt has been made to highlight hard-coded assumptions 
          and segregate these to the headers of classes and the full module. Specific use cases will require modifications.

"""
#from Points import

#Hard-coded assumptions and conversion ratios

default_compression_factors= {
    'default': 0.7,
    'general-siloxane': 0.7,
    }
#Generic dictionaries 
generic_punch_sizes = {
    'default':{
        'diameter_units': 'microns',
        'diameter': 15.875
        },
    'generic_mech_testing':{
        'diameter_units': 'inches',
        'diameter': 5/8
        }
    }
generic_materials = {
    'll50': {},
    'll60': {},
    'SE1700': {},
    'PDMS': {}
    }

# 'print_type'-  'flat'; 'hemi'; 'shaped-special'
generic_volume_dict = {
    'print_type': '',
    'number_of_layers': 0,
    'layer_types':[],
    'strand_sizes': [],
    'layer_arrays': []
    }


#################################################################################################################
###  Lorem  #####################################################################################################
#################################################################################################################

def ideal_volume_guess(structure_dict, 
                       n_structures= 5, 
                       voxel_side_length= 15.875,
                       voxel_resolution_microns= 5,
                       compression_factor= 0.7):
    '''
    Description: Take 

    INPUT:
        'structure_dict'-           lorem
      (optional)
        'n_structures'-             Number of voxels to generate
        'voxel_side_length'-        in mm; total X-Y side length of the full voxel
        'voxel_resolution'-         In microns (um); pixel side length 
        'compression_factor'-       Percentage of 'ideal' thickness 'actual' part will be do to strand-strand convergence pre-curing
    ACTIONS:
    OUTPUTS:
        'return_volume_dict'-       Dictionary with same structure (keys) as 'generic_volume_dict'; represents
                                    results of simulating layer feature 
    '''

    #Initialize variables
    return_volume_dict = {}
    x_array_dim= voxel_side_length *1000 / voxel_resolution_microns   #number of pixels (microns/microns)
    y_array_dim= x_array_dim  #Array is square

    #Check data types and coerce
    if type(compression_factor) == float:
        # check if a fraction
        if compression_factor <=1 and compression_factor >0:
            pass
        else:
            #TODO: add flag for non-fractional compression factors
            pass
    elif type(compression_factor) == int:
        if compression_factor >1 and compression_factor <=100:
            compression_factor = round(compression_factor/100, 5)
        elif compression_factor == 1:
            #Assume no compression and just pass as 1; otherwise compression will be unphysical
            pass

    #Create 



