"""
Copyright 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   20YY-01-26
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Short script to open a digital logbook instance as a Pandas DataFrame.

"""

#Import libraries and functions
import os
import pandas as pd
from tkinter import Tk, filedialog
  # other LayerOver imports
from LayerOver.DataAnalysis.MechanicalData import open_logbook
from LayerOver.DataAnalysis.MechanicalData import split_filename_for_printname_guessing

#String to check for
check_string = "20240813-SJS-01-C.csv"

#Select the file to open
root = Tk()
logbook_filepath = filedialog.askopenfilename()
root.destroy()

#Open the file as a pd.DataFrame
logbook_df = open_logbook(logbook_filepath,
                          target_excel_sheetname = None)

#Check against logbook
printname_guess_dict = split_filename_for_printname_guessing(check_string, logbook_df)
    # Dict keys:
    # 'full_input_name': filename,
    # 'best_guess_printname': '',
    # 'printname_parse_bool': False,
    # 'iteration_marker': '',
    # 'all_potential_matches': [],
    # 'all_match_scores': []

matches = printname_guess_dict['all_matches']
scores = printname_guess_dict['all_match_scores']
best_guess = printname_guess_dict['best_guess_printname']

print(f"The string '{check_string}' has the following matches:")
for match, score in zip(matches, scores):
    print(f"\t {match} \t\t\t {score}")
print(f"\n Best printname guess: {best_guess} \n")
