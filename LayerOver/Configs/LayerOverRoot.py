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

Description: Configure a user root folder for LayerOver settings and data.

TODO:
    -Finalize .layerover root creation
    -Finalize .layerover root check

"""

import os


def create_layerover_root(alt_root = None,
                          working_dir = None):
    '''
    Description: 

    INPUT:
    ACTIONS:
    OUTPUT:
        'generation_dict'       dict; contains filepath and metadata for creation
    '''

    #Initialize filepaths and variables
    root_home = str(Path.home())
    current_workingdir = root_home
    generation_dict = {
        'root_status': 'existing',
        'root_home':'',
        'is_main_root': True,
        'current_working_directory': '',

        }
    
    #Set flags for SEAM repositories
    recipe_check = False
    meta_check = False

    #Apply optional inputs
    if alt_root and os.path.isdir(alt_root):
        root_home = alt_root
    elif os.path.isfile(alt_root):
        # if the supplied alternate fileapath is a file, take the immediate parent directory instead
        if os.path.isdir(os.path.dirname(alt_root)):
            root_home = os.path.dirname(alt_root)
    
    if working_dir and os.path.isdir(working_dir):
        current_workingdir = working_dir
    elif os.path.isfile(working_dir):
        # if the supplied working directory fileapath is a file, take the immediate parent directory instead 
        if os.path.isdir(os.path.dirname(working_dir)):
            current_workingdir = os.path.dirname(working_dir)

    #no matter what, try to find the root LayerOver repo; make one in the root if you can't find it there
    trial_dirs = []
    for root, dirs, files in os.walk(root_home, topdown=False):
        for name in dirs:
            if '.layerover' in name.lower():
                trial_dirs.append(os.path.join(root, name))




        # PRIOR EXAMPLE OF ROOT CREATION FROM SEAM

        # # Check for SEAM root folder and get global
        # root_home = str(Path.home())
        # current_workingdir = root_home
        
        # if alt_root != None:
        #     user_root = alt_root
        #     root_home = self.user_root
        # else:
        #     root_home 
        
        # # .seam root directory and branch
        # if work_dir != None:
        #     self.current_workingdir = work_dir
        
          
        
        #   #alternate directory '.seam' search method; not sure how it's getting to directories outside root          
        # #for item in os.listdir(self.user_root):
        # #    if os.path.isfile(os.path.join(self.user_root, '.seam')):
        # #        trial_dirs.append(os.path.join(self.user_root, '.seam'))

        # if len(trial_dirs) == 1:
            
        #     seam_guess = trial_dirs[0]
            
        #     #Just check that two of the required directorie exist and assume all else is groovy
        #     meta_check = os.path.isdir(os.path.join(seam_guess, "meta_dicts"))
        #     typehash_check = os.path.isdir(os.path.join(seam_guess, "meta_dicts"))
            
        #     if meta_check and typehash_check:
        #         self.seam_root = seam_guess
        #         seam_root = os.path.join(root_home, ".seam") #set the global 'seam_root' as well
        #     else:
        #         self.create_seam_directory(os.path.join(root_home, ".seam"))
        
        # elif len(trial_dirs) == 0:
        #     folder_check = os.path.isdir(os.path.join(root_home, ".seam"))
        #     if not folder_check:
        #         self.create_seam_directory(os.path.join(root_home, ".seam"))
                
        #    #If META or REC filetypes are found, parse into dictionary
    
        #    #Otherwise start a dictionary and save a placeholder (i.e. blank) META and REC file
        #       #Look in SEAM root file 
              
              
        # #set the global 'seam_root' to the
        # seam_root = os.path.join(root_home, ".seam")
        # self.seam_root = seam_root
        
        # self.project_dict = {
        #     '':''
        #     }


def update_structure_data(structure_name, data_directory):

    pass