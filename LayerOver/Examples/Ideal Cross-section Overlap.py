"""
Copyright 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2026-06-06
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Small script to demonstrate cross-section overlap behavior and compression calculation behavior.

Notes:
    - Bottom strand is always the same orientation, with full-thickness at y=0 across all X values.
    - Top strand is offset by 'offset_angle' degrees from the bottom strand, and slice values calculated for each Y-value at that slice's X-value.

"""

import numpy as np
import math
import matplotlib.pyplot as plt

#Direct user-defined variables
strand_diameter = 250                   # in microns
offset_angle = 0
offset_values = [0, 0.1, 0.25, 0.5, 1, 1.5, 2]                # 0<=float<=1; fraction of strand radius to offset; only applies if offset_angle != 0
compression_factor = 0.78

  # other options
overlap_type = 'log_cabin'              # 'log_cabin'- parallel, overlappin strands
strand_thickness_type = 'cylindrical'   # 'cylindrical'
resolution = 1000                       # number of points to simulate across overlap
simulation_range = strand_diameter     # in microns; starting at centerline of bottom strand; should be >= strand_radius to capture full overlap behavior

#Condition variables and functions
if strand_thickness_type == 'cylindrical':
    strand_thickness_func = lambda distance: 2* math.sqrt(strand_radius**2 - distance**2)
strand_radius = strand_diameter/2
  # condition offset values
temp_offsets = []
for offset_value in offset_values:
    if offset_value >= 0 and offset_value <= 2:
        this_value = offset_value*strand_radius
        temp_offsets.append(this_value)
    else:
        temp_offsets.append(offset_value)
offset_values = temp_offsets
  # adjust max height for compression factor
# max_height = strand_diameter*(1-((1-compression_factor)/2))  #half of compression applied to this interface
max_height = strand_diameter*compression_factor 

#Basic parallel offset
for top_offset in offset_values:
    general_layer = np.linspace(0, simulation_range, resolution)
    bottom_layer = np.zeros(len(general_layer))
    top_layer = np.zeros(len(general_layer))
    for idx, x in enumerate(general_layer):
        if x <= strand_radius:
            if strand_thickness_func(x) > 1e-9:
                bottom_layer[idx] = strand_thickness_func(x)
            else:
                bottom_layer[idx] = 0
        top_r = abs(x - top_offset)
        if top_r <= strand_radius:
            top_layer[idx] = strand_thickness_func(top_r)
        else:
            top_layer[idx] = 0

    # modify for max height
    adjusted_top = top_layer.copy()
    adjusted_bottom = bottom_layer.copy()
    adjusted_top[adjusted_top>max_height] = max_height
    adjusted_bottom[adjusted_bottom>max_height] = max_height
    # overlay layers
    overlap_layer = bottom_layer + top_layer
    # adjusted_overlap_diff = ((bottom_layer - adjusted_bottom) + (top_layer- adjusted_top))/2    #divide by 2 to account for overlap of one interface only
    adjusted_overlap_diff = ((bottom_layer - adjusted_bottom) + (top_layer- adjusted_top))
    # create masks
    adjusted_layer = overlap_layer.copy()
    overlap_mask = (top_layer>0) & (bottom_layer>0)                       #raw overlap mask
    # overlap_mask = (top_layer>max_height) & (bottom_layer>max_height)   #overlap of only "high points" that should coallesce
    # difference_layer = np.ones(overlap_layer.shape) * max_height*2
    # difference_layer = overlap_layer - difference_layer
    # difference_layer[difference_layer<0] = 0
    adjusted_layer[overlap_mask] = adjusted_layer[overlap_mask]-adjusted_overlap_diff[overlap_mask]
    report_compress_diff = adjusted_overlap_diff.copy()
    report_compress_diff[~ overlap_mask] = 0

    fig, axs = plt.subplots(2, 1, layout='constrained')
    plt.title(f"Top strand offset by {round(top_offset/strand_radius, 2)}x strand radius")
    axs[0].plot(general_layer, bottom_layer, general_layer, top_layer, general_layer, report_compress_diff)
    #axs[0].set_xlim()
    
    axs[0].set_ylabel('Strand Thickness (microns)')
    axs[0].grid(True)
    axs[0].legend(['Bottom Strand', 'Top Strand', 'Smoosh Adjustment'])

    axs[1].plot (general_layer, overlap_layer, general_layer, adjusted_layer)
    axs[1].set_ylabel('Thickness')
    axs[1].legend(['Thickness sum', 'Adjusted Sum'])
    axs[1].set_xlabel('Distance from bottom-strand centerline (microns)')

    plt.show()

