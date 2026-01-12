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

Description: Small script utilizing LayerOver modules to parse and open mechanical data files.

"""

#Import outside libraries
from tkinter import Tk, filedialog
import os
import matplotlib.pyplot as plt
  
#Imports from LayerOver 
from LayerOver.DataAnalysis import MechanicalData as mech

#%% 
#Select a file(s) to import
root = Tk()
filenames = filedialog.askopenfilenames(filetypes = [('CSV Files', '*.csv')])
root.destroy()

#
for filename in filenames:
    file_output_dict = mech.parse_mech_data_fromcsv(filename)\
    # keys in 'file_output_dict':
    #     'parse_status'                dict
    #         'file_pandas_readable'    bool
    #         'namerow_good'            bool
    #         'namerow_idx'             int
    #     'filename'                    str
    #     'filepath'                    str
    #     'dataframe'                   pandas.DataFrame (hopefully)

    data_df = file_output_dict['dataframe']
    column_names = list(data_df.columns)
    for column in column_names:
        if 'stress' in column.lower():
            stress_col_name = column
        if 'strain' in column.lower():
            strain_col_name = column

    print()
    print("#"*50)
    print(os.path.basename(filename))
    print()
    print(data_df.head(7))

    #Plot all the curves
    # try:
    #     plt.figure(figsize = (10,10))
    #     plt.scatter(data_df[strain_col_name], data_df[stress_col_name])
    #     plt.title(f"{os.path.basename(filename)}")
    #     plt.show()
    # except:
    #     print()
    #     print("Plot failure (hopefully for obvious reasons)")


    replicate_dict = mech.pull_mechanical_replicates(data_df, data_dict = None)
      # pull keys from returned dict
    replicate_parse_success =  replicate_dict['replicate_parse_success']
    number_of_replicates = replicate_dict['number_of_replicates']
    last_strain_peak_index = replicate_dict['last_strain_peak_index']
    last_strain_valley_index = replicate_dict['last_strain_valley_index']
    peak_to_valley_index_diff = replicate_dict['peak_to_valley_index_diff']
      # each value in "replicate_data_dict" is a pandas.DataFrame (hopefully)
      # replicate numbering starts at 1
    replicate_data_dict = peak_to_valley_index_diff = replicate_dict['replicate_data']
    try:
        last_df = replicate_data_dict[number_of_replicates]

        print()
        print("#"*50)
        print(os.path.basename(filename))
        print()
        if replicate_parse_success:
            print(f"Number of replicates parsed: {number_of_replicates}")
        else:
            print("Failure to pull mechanical data replicates.")
        print()
    except KeyError:
        print()
        print("Failure to parse file (no replicate mechanical data found)")

    try:
        plt.figure(figsize = (10,10))
        plt.scatter(last_df['strain_data_loading'], last_df['stress_data_loading'], color = 'r')
        plt.scatter(last_df['strain_data_unloading'], last_df['stress_data_unloading'], color = 'g')
        plt.title(f"{os.path.basename(filename)}")
        plt.show()
    except:
        print()
        print("Plot failure (hopefully for obvious reasons)")