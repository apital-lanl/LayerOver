"""
Copyright 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   20YY-MM-DD
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Short script to open a digital logbook instance as a Pandas DataFrame.

"""

#Import libraries and functions
import os
import pandas as pd
from tkinter import Tk, filedialog
  # other LayerOver imports
from LayerOver.Analysis.MechanicalData import open_logbook
from LayerOver.Analysis.MechanicalData import split_filename_for_printname_guessing

#Select the file to open
root = Tk()
logbook_filepath = filedialog.askopenfilename()
root.destroy()

#Open the file as a pd.DataFrame
logbook_df = open_logbook(logbook_filepath,
                          target_excel_sheetname = None)

#Get a glance at the logbook
print(logbook_df.head(50))
print()

#Check column names
print(f"Columns: ")
for column_name in list(logbook_df.columns):
    print(f"\t {column_name}")

#Save DataFrame to CSV with a similar name
savebook_filename = os.path.basename(logbook_filepath).replace('AutomatedAnalysis', 'AutomatedResults')
savebook_filepath = os.path.join(os.path.dirname(logbook_filepath), savebook_filename)
with open(savebook_filepath, 'w') as file:
    logbook_df.to_csv(file, index= False, lineterminator='\n')


