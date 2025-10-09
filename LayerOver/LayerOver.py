
# -*- coding: utf-8 -*-
"""
© 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2024-09-12
Modified:  2025-07-07
Version:   1.2.0 

@author: Aaron Pital (Los Alamos National Lab)

Description: Combined library for gcode, STL, voxelization, and related operations for Direct Ink Write (DIW) work.

"""
version = '1.2.0'
last_modified_date = '2025-07-07'

  #System and built-ins
import os
import math
from tkinter import Tk, filedialog
import traceback

  #Visualizaiton
import matplotlib.pyplot as plt

  #Data Handling
import numpy as np

  #Utilities
from tqdm.auto import tqdm

  #LayerOver imports
from Part import Gcode
from Points import Point_Clod
import Voxel
import Opacity


# Define package-variables
color_lists = {
'Rainbow':['salmon','darkorange','gold','limegreen',\
           'mediumturquoise','slateblue','mediumorchid'],
'rRainbow': ['salmon','darkorange','gold','limegreen',\
           'mediumturquoise','slateblue','mediumorchid'].reverse(),
}
    
#################################################################################################################

        



