"""
Copyright 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-01-20
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: 

"""

from LayerOver.Analysis import MechanicalData as mech
import os
import pandas as pd
pd.set_option('display.max_columns', None)  #Make displayed columns full-width
from tkinter import Tk, filedialog

#Select the directory
root = Tk()
mech_directory = filedialog.askdirectory(title = "Select the directory with stress-strain data")
root.destroy()

#Do the stuff
mech_data_df = mech.process_directory_for_mech_files(directory = mech_directory, 
                                                     show_each_file_results = False,
                                                     save_last_replicate_graph = True)

# #Show the stuff
# print(f"Quick-look at mech data from {mech_directory}")
# print(f"{mech_data_df.head(25)}")
# print()

# #Save the stuff
# #Commented out because this is done in 'process_directory_...' already
# save_dir = os.path.join(mech_directory, 'Extracted Mechanical Summary Data')
# if os.path.isdir(save_dir):
#     save_name = "Extracted Mechanical Data.csv"
#     save_path = os.path.join(save_dir, save_name)
#     mech_data_df.to_csv(save_path)
# else:
#     os.makedirs(save_dir, exist_ok=False)
#     save_name = "Extracted Mechanical Data.csv"
#     save_path = os.path.join(save_dir, save_name)
#     mech_data_df.to_csv(save_path)