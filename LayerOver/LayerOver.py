# -*- coding: utf-8 -*-
"""
Created:   2024-09-12
Modified:  2025-04-17
Version:   1.1.0 

@author: Aaron Pital (Los Alamos National Lab)

Description: Combined library for gcode, STL, voxelization, and related operations for Direct Ink Write (DIW) work.

"""
version = '1.1.0'
last_modified_date = '2025-04-17'

  #System and built-ins
import os
import sys
import json
import math
import csv
from tkinter import Tk, filedialog

  #Visualizaiton
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

  #Data Handling
import numpy as np
import pandas as pd

  #Scientifiic algorithm packages
import scipy.interpolate
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

  #Get random biz
from random import choice
import random

  #Utilities
from tqdm.auto import tqdm


# Define package-variables
color_lists = {
'Rainbow':['salmon','darkorange','gold','limegreen',\
           'mediumturquoise','slateblue','mediumorchid'],
'rRainbow': ['salmon','darkorange','gold','limegreen',\
           'mediumturquoise','slateblue','mediumorchid'].reverse(),
}
    

###############################################################################################################################
    
        
        
###############################################################################################################################       
    

    
    
###############################################################################################################################
    
    

###############################################################################################################################

#%%

