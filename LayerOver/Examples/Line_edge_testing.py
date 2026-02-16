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
import math
import matplotlib.pyplot as plt
from LayerOver.Core.Points import get_radial_neighbors_array_data

interior_points = (20,45)
angle = 90
line_radius = 10
array_dim = (150, 150)
length = None
line_type = 'simple'
strand_diameter = 50

interior_point_list = [
    [20, 30],
    [2,6],
    [17,149],
    [149, 30],
    ]
# angle_list = [
#     0,1,40,90,135,179,180]
angle_list = [115, 135, 145, 170]
    
"""
Description:
    Interpret line points as a strand of radius 'line_radius' and draw in an array
INPUTS:
    'interior_point'    array, tuple, or list of iterables; interpret as tuple of array indices (Y,X)
    'angle'             float or int; angle between 'interior_point' and array horizontal axis
                        if 'length' is None or 0, assume line is drawn across entire array with center at 'interior_point'
                        otherwise, assume line is drawn in only one direction from 'interior_point' at 'angle'
    'array_dim'         int, tuple, or list; dimensions of 2D array to draw line on
    (optional)
    'line_radius'       float or int; stand radius in number of pixels (i.e. array indices distance)
                        NOTE: thickness will be 2*'line_radius' at center of line, so by default 3 lines will be drawn on average; centerline and parallel radial lines 1 pixel over, one on each side
    'length'            float or int; if NONE, assume line fills the screen edge-to-edge; othewise, assume 'interior_point' is the starting point
                        NOTE: only option for drawing a line in both directions is if 'length'= None. Otherwise 'angle' determines a 2-D line direction.
    'line_type'         str; 'simple'- 2D line from point-to-point
                        'bezier'- bezier curve; control points will be interpolated
ACTIONS:
    -lorem
OUTPUTS:
    'control_point_dict'    dict; specifics of line end points and the 
"""
for interior_points in interior_point_list:
    for angle in angle_list:

        #Initialize variables
        control_point_dict = {
            'control_points': {},
            'interior_only': True,
            'drawn_array': None
            }
        point_dict = {
            'coordinates': [],
            'radius_values': [],

            }
        start_points = []
        end_points = []
    
        #Condition input coordinates
        if type(interior_points) == tuple:
            interior_points = np.array([[interior_points[0], interior_points[1]]])
        if type(interior_points) == list:
            #TODO: check that list elements are actually (Y,X) formatted
            if type(interior_points[0])== list:
                interior_points = np.array(interior_points)
            elif (type(interior_points[0]==int) or (type(interior_points[0]==float))):
                interior_points = np.array([interior_points])
            # coerce to int because these are array indices
        interior_points = interior_points.astype(np.int32)
    
        #Make sure angle is appropriate
        if angle >360:
            coerced_angle = angle
            while coerced_angle >360:
                coerced_angle = coerced_angle%360
            angle= coerced_angle

        #Condition array
        if (type(array_dim) == int) or (type(array_dim) == float):
            array_x_dim = int(array_dim)
            array_y_dim = int (array_dim)
        elif type(array_dim) == tuple:
            if len(array_dim) == 2:
                array_x_dim = array_dim[1]
                array_y_dim = array_dim[0]
            else:
                print()
                print("Bad input values for 'array_dim'; check values")
                print()
        elif type(array_dim) == list:
            if len(array_dim) == 2:
                array_x_dim = array_dim[1]
                array_y_dim = array_dim[0]
            else:
                print()
                print("Bad input values for 'array_dim'; check values")
                print()
        return_array = np.zeros((array_y_dim, array_x_dim))

        #Set conditions for drawing lines
            # set edge-finding conditions
        if (not length) and (line_type == 'simple'):
            #In this case, assume a line should be drawn across the entire screen
            control_point_dict['interior_only'] = False
            #Find closest array edges
            for idx, point in enumerate(interior_points):
                this_y = point[0]
                this_x = point[1]
                #Get upper line angle
                if angle >180:
                    upper_angle = angle%180
                else:
                    upper_angle = angle
                #Get slope of line
                slope = math.tan(math.radians(upper_angle))
                #Get max_y at max_y and vice versa
                  # X is easy because it's always > to the right; 
                y_at_right = round(((array_x_dim-this_x) * slope *-1) + this_y)  #Flip sign of slope for delta-y
                y_at_left =  round(((this_x) * slope) + this_y)
                  # Slope changes Y behaviour (i.e. how much 'Y' is left for each side), so both cases have to be accounted for
                if slope > 0:
                    if abs(slope) > 1e-7:
                        x_at_right = round(this_x + abs((this_y)/slope) )
                        x_at_left = round(this_x - abs((array_y_dim- this_y)/slope))
                    #If slope is basically zero, just assume horizontal to avoid division by ~zero overflow error
                    else:
                        x_at_left = 0
                        x_at_right = (array_x_dim - 1)
                elif slope < 0:
                    if abs(slope) > 1e-7:
                        x_at_right = round(this_x + abs((array_y_dim- this_y)/slope) )
                        x_at_left = round(this_x - abs((this_y)/slope))
                    #If slope is basically zero, just assume horizontal to avoid division by ~zero overflow error
                    else:
                        x_at_left = 0
                        x_at_right = (array_x_dim - 1)
                else:
                    x_at_left = 0
                    x_at_right = (array_x_dim - 1)
                
                #Check to see if edge indices are within array bounds
                y_right_bool = bool((y_at_right >=0) and (y_at_right<=(array_y_dim-1)))
                x_right_bool = bool((x_at_right >=0) and (x_at_right<=(array_x_dim-1)))
                y_left_bool = bool((y_at_left >=0) and (y_at_left<=(array_y_dim-1)))
                x_left_bool = bool((x_at_left >=0) and (x_at_left<=(array_x_dim-1)))

                #Find appropriate array edge coordinates
                  # for right side of point
                if y_right_bool and x_right_bool:
                    #If both edges have valid indices, find the closest one
                    right_y_guess = [y_at_right ,(array_x_dim-1)]
                    if slope>0:
                        right_x_guess = [0, x_at_right]
                    elif slope<0: 
                        right_x_guess = [(array_y_dim-1), x_at_right]

                    #Distance to top/bottom and right edges
                    x_edge_distance = math.sqrt((this_x-right_x_guess[1])**2 + 
                                           (this_y-right_x_guess[0])**2)
                    y_edge_distance = math.sqrt((this_x-right_y_guess[1])**2 + 
                                           (this_y-right_y_guess[0])**2)

                    if x_edge_distance<y_edge_distance:
                        right_edge_point = right_x_guess
                    elif x_edge_distance > y_edge_distance:
                        right_edge_point = right_y_guess
                elif y_right_bool:
                    right_edge_point = [y_at_right, (array_x_dim-1)]
                elif x_right_bool:
                    if slope > 0:
                        right_edge_point = [0, x_at_right]
                    else:
                        right_edge_point = [(array_y_dim-1), x_at_right]
                  # for the left side of the point
                if y_left_bool and x_left_bool:
                    #If both edges have valid indices, find the closest one
                    left_y_guess = [y_at_left ,(array_x_dim-1)]
                    if slope<0:
                        left_x_guess = [0, x_at_left]
                    elif slope>0: 
                        left_x_guess = [(array_y_dim-1), x_at_left]

                    #Distance to top/bottom and right edges
                    x_edge_distance = math.sqrt((this_x-left_x_guess[1])**2 + 
                                           (this_y-left_x_guess[0])**2)
                    y_edge_distance = math.sqrt((this_x-left_y_guess[1])**2 + 
                                           (this_y-left_y_guess[0])**2)

                    if x_edge_distance < y_edge_distance:
                        left_edge_point = left_x_guess
                    elif x_edge_distance > y_edge_distance:
                        left_edge_point = left_y_guess
                elif y_left_bool:
                    left_edge_point = [y_at_left, 0]
                elif x_left_bool:
                    if slope >0:
                        left_edge_point = [(array_y_dim-1), x_at_left]
                    else:
                        left_edge_point = [0, x_at_left]
        

        print(f"For {interior_points[0]} at angle of {angle}:")
        print(f"\t y_at_right: {y_at_right} \t\t {y_right_bool}")
        print(f"\t x_at_right: {x_at_right} \t\t {x_right_bool}")
        print(f"\t y_at_left: {y_at_left} \t\t {y_left_bool}")
        print(f"\t x_at_left: {x_at_left} \t\t {x_left_bool}")
        print()
        print(f"Left-edge guess: {left_edge_point}")
        print(f"Right-edge guess: {right_edge_point}")
        print('_'*50)
        print()

        ##Plot the line segment 
        # plt.scatter([left_edge_point[1], right_edge_point[1]],[left_edge_point[0], right_edge_point[0]], marker = 'x', color = 'r')
        # plt.scatter([this_x], [this_y], marker='o', color = 'k')
        # plt.plot([left_edge_point[1], right_edge_point[1]],[left_edge_point[0], right_edge_point[0]], color = 'r')
        # plt.title(f"Point {interior_points[0]} at angle of {angle}")
        # plt.xlim(0, array_x_dim)
        # plt.ylim(array_y_dim, 0)
        # plt.show()

        #Test radial line getter for plotted line
        radial_dict = get_radial_neighbors_array_data(left_edge_point,
                                                       right_edge_point,
                                                       angle,
                                                       strand_diameter,
                                                       array_dim,
                                                       pix_to_um_conv = 1,
                                                       thickness_fcn = 'cylinder',
                                                       print_intermediate_steps = True)
        start_points = radial_dict['starting_points']
        end_points = radial_dict['ending_points']
        distances = radial_dict['thicnkess']

        print()
        plt.scatter([left_edge_point[1], right_edge_point[1]],[left_edge_point[0], right_edge_point[0]], marker = 'x', color = 'r')
        plt.scatter([this_x], [this_y], marker='o', color = 'k')
        plt.plot([left_edge_point[1], right_edge_point[1]],[left_edge_point[0], right_edge_point[0]], color = 'r')
        for start, stop, distance in zip(start_points, end_points, distances):
        #     print(f"\t {start} \t {stop} \t {distance}")
            plt.plot([start[1], stop[1]],[start[0], stop[0]], color = 'gray', linewidth = 1, alpha = 0.8)
        print('#'*50)

        plt.title(f"Point {interior_points[0]} at angle of {angle}")
        plt.xlim(0, array_x_dim)
        plt.ylim(array_y_dim, 0)
        plt.show()