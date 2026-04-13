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

Description: 

"""

import numpy as np
from LayerOver.Core.Points import draw_2D_strand_line_bythickness

#Main user defined parameters
show_intermediate_step_plots = True
print_intermediate_steps = True
strand_diam = 150           #in microns, converted to pixels below
array_dim = (5000, 5000)
array_side_length_mm = 15   #in mm

#Secondary user defined parameters
# interior_points = (20,45)
# angle = 90
voxel_resolution_microns = round((array_side_length_mm*1000)/max(array_dim), 5)
strand_diam_pix = strand_diam/voxel_resolution_microns
length = None
line_type = 'simple'
strand_diameter = 50
  # testing parameters
interior_point_list = [
    [200, 300],
    [20,60],
    [170,1490],
    [1490, 300],
    ]
# angle_list = [
#     0,1,40,90,135,179,180]
angle_list = [20, 45, 115, 135, 170]


for angle in angle_list:
    for interior_points in interior_point_list:
        control_point_dict = draw_2D_strand_line_bythickness(interior_points,
                                                    angle,
                                                    strand_diam_pix,
                                                    array_dim,
                                                    um_to_pix_conversion = voxel_resolution_microns,
                                                    length = None,
                                                    line_type = 'simple',
                                                    thickness_fcn = 'cylinder',
                                                    show_points = False,
                                                    show_array_iterations = show_intermediate_step_plots,
                                                    show_final_array = True, 
                                                    report_line_drawing_iterations = print_intermediate_steps)