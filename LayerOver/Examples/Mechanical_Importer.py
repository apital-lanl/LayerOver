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

Description: Small script utilizing LayerOver modules to parse and open mechanical data files

"""

#Import outside libraries
from tkinter import Tk, filedialog
import os
import matplotlib.pyplot as plt
  
#Imports from LayerOver 
from LayerOver.DataAnalysis import MechanicalData as mech


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


    last_df = mech.pull_last_mechanical_replicate(data_df, data_dict = None)

    print()
    print("#"*50)
    print(os.path.basename(filename))
    print()
    print(last_df.head(7))

    try:
        plt.figure(figsize = (10,10))
        plt.scatter(data_df[strain_col_name], data_df[stress_col_name])
        plt.title(f"{os.path.basename(filename)}")
        plt.show()
    except:
        print()
        print("Plot failure (hopefully for obvious reasons)")