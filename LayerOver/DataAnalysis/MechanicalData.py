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

Description: Module for pulling and analyzing mechanical data from Excel and CSV files. 

"""

#Import libraries
import os
from tkinter import filedialog, Tk

  # import other LayerOver modules


#Define variables

######################################################################################################################################
#####                #################################################################################################################
######################################################################################################################################

def walk_directory_for_mech_files(directory=None):
    if not directory:
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()

    elif os.path.isdir(directory):
        #Walk the directory


    else:
        print("Bad directory; please select an appropriate directory with mechanical data.")
        root = Tk()
        directory = filedialog.askdirectory(title="Select directory with all the mechanical data.")
        root.destroy()