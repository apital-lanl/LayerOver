
# -*- coding: utf-8 -*-
"""
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

        



