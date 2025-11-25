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


class CompressionFactors(self):

    #Hard-coded assumptions and conversion ratios
    generic_punch_sizes = {
        'generic_mech_testing':{
            'diameter_units': 'inches',
            'diameter': 5/8
            }
        }

    materials = {
        'll50': {},
        'll60': {},
        'SE1700': {},
        'PDMS': {}
        }