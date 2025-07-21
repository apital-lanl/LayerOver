# -*- coding: utf-8 -*-
"""
Created:   2024-04-10
Modified:  2024-05-12
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: Object class for POV renders, POV plots, and geometry-based measures of a collection of XYZ cartesian geometries.

"""
version = '0.0.0'
last_modified_date = '2000-00-00'

  #System and built-ins
import math
import os
import sys
from random import choice
from tkinter import Tk, filedialog
import traceback

  #Data Handling
import numpy as np

  #Visualizaiton
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits import mplot3d

  #STL and meshing
from skimage.measure import marching_cubes
from stl.mesh import Mesh

  #Utilities
from tqdm.auto import tqdm

  #Other LayerOver modules
import Points
from Part import Gcode


def get_random_point_summaries(gcode_part, n_retrieval_points = 5, \
                               excluded_point_indcs = [], show_iterations = False):

    ''' v0.2.2  created:2024-05-08  modified:2024-05-12 
    
    INPUT:   'gcode_part'- a part generated from 'gcode' class; inclueds config.JSON and .PGM layer files
        (optional)
             'n_retrieval_points' - number of points to evaluate and return
             
    ACTION:  Pull fibonacci point triangles from part and calculate center, normal, and other stats 
    OUTPUT:  Dictionary with each point's summary data; dictionary keys are just integer #'s starting at 0 (i.e. index of point_
    '''
    
    #v0.2.0 version; takes known fibonacci points (that seem like they work) and makes a quick 'nearest neighbor' triangle to other fib. points
      #Get stats about part and initialize variables
    point_n_total = gcode_part.substrate_dict['fib_points'].shape[0]
    eval_point_dict = {}
    
    #Get 'n' triangles from the surface randomly
    rand_point_indcs = []
    rand_triang_coord = []
    for i in range(0, n_retrieval_points):
        this_choice = choice(range(0, point_n_total))
        this_point = gcode_part.substrate_dict['fib_points'][this_choice]
          #add to list if it's not already there (don't want duplicates)
        these_coords = [this_point]
        if (this_choice not in rand_point_indcs) and (this_choice not in excluded_point_indcs):
            rand_point_indcs.append(this_choice)
            trial_points=[]
            for idx, trial_point in enumerate(gcode_part.substrate_dict['fib_points']):
                this_distance = math.sqrt((trial_point[0]-this_point[0])**2 + (trial_point[1]-this_point[1])**2 + (trial_point[2]-this_point[2])**2)
                if this_distance != 0:
                    trial_points.append((idx, this_distance))
            #Sort coordinates by distance and add two nearest points
            sorted_trial_points = sorted(trial_points, key=lambda x: x[1])
              #smallest values sorted first; pick the smallest two
            second_point = gcode_part.substrate_dict['fib_points'][sorted_trial_points[0][0]]
            these_coords.append(second_point)
              #calculate angle betwen p0->p1 and p0->p2; if about 180 (or 0), use a new p2
            third_point = gcode_part.substrate_dict['fib_points'][sorted_trial_points[1][0]]
            vector_1 = second_point-this_point
            vector_2 = third_point-this_point
            this_angle = math.acos( ((vector_1[0]*vector_2[0]+vector_1[1]*vector_2[1]+vector_1[2]*vector_2[2])/  \
                        (math.sqrt(vector_1[0]**2+vector_1[1]**2+vector_1[2]**2) * \
                         math.sqrt(vector_2[0]**2+vector_2[1]**2+vector_2[2]**2))))
            if this_angle >2:
                #try another point; shouldn't be possible to be co-linearish at this point
                third_point = gcode_part.substrate_dict['fib_points'][sorted_trial_points[2][0]]
                vector_1 = second_point-this_point
                vector_2 = third_point-this_point
                this_angle = math.acos( ((vector_1[0]*vector_2[0]+vector_1[1]*vector_2[1]+vector_1[2]*vector_2[2])/  \
                        (math.sqrt(vector_1[0]**2+vector_1[1]**2+vector_1[2]**2) * \
                         math.sqrt(vector_2[0]**2+vector_2[1]**2+vector_2[2]**2))))
                these_coords.append(third_point)
            else:
                these_coords.append(third_point)
              #Add triangle coordinates to list
            rand_triang_coord.append(these_coords)
    
    ''' v0.1 version; uses convex hull function that I don't trust fully
    #Get stats about part and initialize variables
    triangle_n_total = gcode_part.substrate_dict['fib_surface_triangles'].shape[0]
    eval_point_dict = {}
    
    #Get 'n' triangles from the surface randomly
    rand_triang_indcs = []
    rand_triang_points = []
    rand_triang_coord = []
    for i in range(0, n_retrieval_points):
        this_choice = choice(range(0, triangle_n_total))
        these_indcs = gcode_part.substrate_dict['fib_surface_triangles'][this_choice]
          #add to list if it's not already there (don't want duplicates)
        these_coords = []
        if (this_choice not in rand_triang_indcs) and (this_choice not in excluded_point_indcs):
            rand_triang_indcs.append(this_choice)
            rand_triang_points.append(these_indcs)
            for idx in these_indcs:
                these_coords.append(gcode_part.substrate_dict['fib_points'][idx])
            rand_triang_coord.append(these_coords)
    '''
    
    #Run through each set of triangles and get normal, center, equation of plane, etc. for them
    point_dicts = {}
    for idx, point_coords in enumerate(rand_triang_coord):
        this_point_dict = Points.get_3point_normal(point_coords, show_plot = show_iterations)
        point_dicts.update({idx:this_point_dict})

    return point_dicts


def get_layer_points(gcode_part, voxel_dict, boundary_buffer_multiplier = 0.1):
    
    ''' v0.2.0  created:2024-05-09  modified:2024-06-03

    Take a part with layers and a dictionary of voxel properties and return a dictionary of layer-points within
      the voxel region and some properties.
    
        INPUT:   'gcode_part'- gcode class instance with layers added (happens at initialization if 'Layer -' files are selected)
                 'voxel_dict'- dictionary with voxel bounds, etc.
        ACTION:  cull points to within voxel bounds; easy because they're square at this point (that's a joke)
        OUTPUT:  updated 'voxel_dict' with layers and in-bound coordinates added as lists;
                    NOTE: lists have inhomogenous dimensions, so naive parsing is kind of a pain.
    
    '''

    #Initialize new entry for layer coordinates
    voxel_dict.update( {'layer_coordinates':{}})
    
    #TODO: make sure the output of the 'voxel_bounding_box' is ordered somehow
    bounds_array = []
    for key in list(voxel_dict['voxel_bounding_box']):
        bounds_array.append(voxel_dict['voxel_bounding_box'][key])
    xmin, xmax, ymin, ymax, zmin, zmax = bounds_array
    xmin -= abs(xmax-xmin)*boundary_buffer_multiplier
    xmax += abs(xmax-xmin)*boundary_buffer_multiplier
    ymin -= abs(ymax-ymin)*boundary_buffer_multiplier
    ymax += abs(ymax-ymin)*boundary_buffer_multiplier
    zmin -= abs(zmax-zmin)*boundary_buffer_multiplier
    zmax += abs(zmax-zmin)*boundary_buffer_multiplier

    #run through each layer's coordinates and cull to dictionary list
    for layer_name in list(gcode_part.layers.keys()):
          #initialize dictionary entry for this layer
        voxel_dict['layer_coordinates'].update( {layer_name:{}} )
        
        this_layer_dict = gcode_part.layers[layer_name]
        these_coordinates = this_layer_dict['coordinates']
        #save as lists for appending to dictionary
        in_bound_indcs = []
        in_bound_coordinates = []
        seq_coordinates = []
        last_idx_hit = -1  #placeholder for denoting sequential coordinates
        for idx, coordinates in enumerate(these_coordinates):
            x, y, z = coordinates
            
            if (x>xmin) and (x<xmax) and (y>ymin) and (y<ymax) and (z>zmin) and (z<zmax):
                in_bound_indcs.append(idx)
                #want to keep sequential coordinates nested
                if (idx-last_idx_hit) > 1:
                    in_bound_coordinates.append(seq_coordinates)  #add the previous list
                    seq_coordinates = [coordinates]  #start a new list with this index's coordinates
                else:
                    seq_coordinates.append(coordinates)
                    
                last_idx_hit = idx
        in_bound_coordinates.append(seq_coordinates)  #append whatever is left in the queue

        #append lists to dictionary under 'layer_name' for later pulling from gcode
        voxel_dict['layer_coordinates'][layer_name].update( {'indices': in_bound_indcs} )
        voxel_dict['layer_coordinates'][layer_name].update( {'coordinates': in_bound_coordinates} )
            
    return voxel_dict


def get_voxel(center, vector, radius= 3, voxel_depth= 6, n_pixels=100):
    
    ''' v0.2.4  created:2024-05-09  modified:20245-04-17
    v0.2.3- modified bounding points to be actual 'pixel' array rather than radial points
    v0.2.4- added 'top' or 'bottom' voxel face check
    
    Given a point and a vector, generate two circles at two points along with 'center' in middle.
        Generate a set of min:max bounds defining rectangular voxel for that cylindrical region.
        Will be bounding box for region, not true voxel region.
        INPUT:   'vector'- should be normal vector to plane
        ACTION:  lorem
        OUTPUT:  lorem
        
    Changes:
        2025-04-17- changed 'n_circ_points' to reflect change in function. Used to be (360/n_circ_points)= "num of points returned"; now, n_circ_points is simply the "num of points returned"
    '''

    #Norm vector just in case it's not already
    vector /= np.linalg.norm(vector)
    
    #Get two bounding points for cylindrical region (center of each end's face); 'top'/'bottom' are arbitrary
    top_center = center + (vector*(voxel_depth/2))
    bottom_center = center - (vector*(voxel_depth/2))
      #Generate points for circumerence of cirlce
      #   will be 360/'n_circ_points' number of points
    top_radial_points = Points.generate_radial_points(vector, top_center, radius, n_circ_points = 6)
    bottom_radial_points = Points.generate_radial_points(vector, bottom_center, radius, n_circ_points = 6)

    #Get corners of visualization square
    pixel_array_radius = math.sqrt(2*radius**2)
    top_pixel_array_corners = Points.generate_radial_points(vector, top_center, pixel_array_radius, n_circ_points = 4)
    bottom_pixel_array_corners = Points.generate_radial_points(vector, bottom_center, pixel_array_radius, n_circ_points = 4)
    
    #Draw 'top' pixel array by coordinates
    pa = top_pixel_array_corners[0]
    pb = top_pixel_array_corners[1]
    pc = top_pixel_array_corners[2]

    # Generate point-wise grid based on plane equation; 
    pixel_coordinate_array = np.zeros((3, n_pixels**2))
    idx = 0
    for u in tqdm(range(0, n_pixels)):
        for v in range(0,n_pixels):
            this_point = np.multiply((u/n_pixels), pa) + \
                         np.multiply((v/n_pixels), pb) + \
                         np.multiply(((n_pixels-u-v)/n_pixels), pc)
            x,y,z = this_point
            pixel_coordinate_array[0,idx] = x
            pixel_coordinate_array[1,idx] = y
            pixel_coordinate_array[2,idx] = z
            idx += 1
    temp_top_pixel_coordinate_array = pixel_coordinate_array.transpose()
    
    #Draw 'bottom' pixel array by coordinates
    pa = bottom_pixel_array_corners[0]
    pb = bottom_pixel_array_corners[1]
    pc = bottom_pixel_array_corners[2]

    # Generate point-wise grid based on plane equation; 
    pixel_coordinate_array = np.zeros((3, n_pixels**2))
    idx = 0
    for u in tqdm(range(0, n_pixels)):
        for v in range(0,n_pixels):
            this_point = np.multiply((u/n_pixels), pa) + \
                         np.multiply((v/n_pixels), pb) + \
                         np.multiply(((n_pixels-u-v)/n_pixels), pc)
            x,y,z = this_point
            pixel_coordinate_array[0,idx] = x
            pixel_coordinate_array[1,idx] = y
            pixel_coordinate_array[2,idx] = z
            idx += 1
    temp_bottom_pixel_coordinate_array = pixel_coordinate_array.transpose()

    #Check which array is 'outside' (further away from origin) and make that 'top_array'
      #pick first points and measure distance from origin
    ptopx, ptopy, ptopz = temp_top_pixel_coordinate_array[0]
    pbotx, pboty, pbotz = temp_bottom_pixel_coordinate_array[0]
    top_distance = math.sqrt(ptopx**2 + ptopy**2 + ptopz**2)
    bot_distance = math.sqrt(pbotx**2 + pboty**2 + pbotz**2)
    if bot_distance > top_distance:
        bottom_pixel_coordinate_array = temp_top_pixel_coordinate_array
        top_pixel_coordinate_array = temp_bottom_pixel_coordinate_array
        print("Flipped: ")
        print(f"Top: {temp_bottom_pixel_coordinate_array[0]}  {bot_distance}")
        print(f"Bot: {temp_top_pixel_coordinate_array[0]}  {top_distance}")
        print()
    else:
        bottom_pixel_coordinate_array = temp_bottom_pixel_coordinate_array
        top_pixel_coordinate_array = temp_top_pixel_coordinate_array
        print("Not Flipped: ")
        print(f"Top: {temp_top_pixel_coordinate_array[0]}  {top_distance}")
        print(f"Bot: {temp_bottom_pixel_coordinate_array[0]}  {bot_distance}")
        print()
        
    #Whip through circumferential points and get min:max bounding box boundaries 
    x_min = 10000
    x_max = -10000
    y_min = 10000
    y_max = -10000
    z_min = 10000
    z_max = -10000
    
    #used to be 'top_radial_points'
    for point in top_pixel_coordinate_array:
        if point[0] < x_min:
            x_min = point[0]
        if point[0] > x_max:
            x_max = point[0]
        if point[1] < y_min:
            y_min = point[1]
        if point[1] > y_max:
            y_max = point[1]
        if point[2] < z_min:
            z_min = point[2]
        if point[2] > z_max:
            z_max = point[2]
    for point in bottom_pixel_coordinate_array:
        if point[0] < x_min:
            x_min = point[0]
        if point[0] > x_max:
            x_max = point[0]
        if point[1] < y_min:
            y_min = point[1]
        if point[1] > y_max:
            y_max = point[1]
        if point[2] < z_min:
            z_min = point[2]
        if point[2] > z_max:
            z_max = point[2]

    #Make a return_dict
    voxel_dict = {
        'top_cylinder_circum_points': top_radial_points,
        'bottom_cylinder_circum_points': bottom_radial_points,
        'top_pixel_array_corners':  top_pixel_array_corners,
        'bottom_pixel_array_corners':  bottom_pixel_array_corners,
        'top_pixel_coordinate_array': top_pixel_coordinate_array,
        'bottom_pixel_coordinate_array': bottom_pixel_coordinate_array,
        #TODO: add plane equations
        #point_clod.get_3point_normal(point_coords)
        'voxel_bounding_box': {
            'xmin': x_min,
            'xmax': x_max,
            'ymin': y_min,
            'ymax': y_max,
            'zmin': z_min,
            'zmax': z_max},
        'center': center,
        'vector': vector,
        'voxel_vertices': {
             1: [x_min, y_min, z_min],
             2: [x_max, y_max, z_min],
             3: [x_min, y_max, z_min],
             4: [x_max, y_min, z_min],
             5: [x_min, y_min, z_max],
             6: [x_max, y_max, z_max],
             7: [x_min, y_max, z_max],
             8: [x_max, y_min, z_max]
                },
        'voxel_faces': {
            'bottom':[1,3,2,4],
            'left':[1,3,7,5],
            'top':[5,8,6,7],
            'front':[1,4,8,5],
            'right':[4,2,6,8],
            'back':[2,6,7,3]
            }
                }

    return voxel_dict


def voxel_points_visualize(voxel_dict, 
                       gif = False,
                       gif_frames = 4,
                       gif_integrals = 50,
                       azim_total_degrees = 360,
                       azim_step = 0,
                       save_plot = True):

    '''
    INPUT:
        gif_frames          Number of total frames in the gif
        gif_integrals       Delay between frames of gif (in msec)
        azim_total_degrees  Number of degrees to spin the graph during the gif
        azim_step           Number of degrees to step spinning graph in each frame of gif
    OUTPU:
        Returns nothing
        Opens new window with plot, plots inline, a 
    '''
            
    plot_first_layer = True
    plot_last_layer = True
    
    if azim_step == 0:
        azim_step = azim_total_degrees/gif_frames
            
    top_pixel_coordinate_array = voxel_dict['top_pixel_coordinate_array']
    bottom_pixel_coordinate_array = voxel_dict['bottom_pixel_coordinate_array']
    layer_list = voxel_dict['layer_names']
    voxel_name = voxel_dict['voxel_name']
    save_name = voxel_dict['save_filepath']
            
    bottom_layer_name = layer_list[0]
    top_layer_name = layer_list[-1]
            
    first_points = np.array([0])
    for c_idx, contig_list in enumerate(voxel_dict['layer_coordinates'][bottom_layer_name]['coordinates']):
        for p_idx, point in enumerate(contig_list):
            this_point = point
            #ignore 0-dim arrays
            if this_point.shape[0]<3:
                pass
            elif len(first_points.shape)<2:
                first_points = np.array([this_point])
            else:
                first_points = np.append(first_points, [this_point], axis =0)
    first_points = np.transpose(first_points)
            
    last_points = np.array([0])
    for c_idx, contig_list in enumerate(voxel_dict['layer_coordinates'][top_layer_name]['coordinates']):
        for p_idx, point in enumerate(contig_list):
            this_point = point
            #ignore 0-dim arrays
            if this_point.shape[0]<3:
                pass
            elif len(last_points.shape)<2:
                last_points = np.array([this_point])
            else:
                last_points = np.append(last_points, [this_point], axis =0)
    last_points = np.transpose(last_points)
            
        #split x/y/z for plotting
    top_xs = top_pixel_coordinate_array[::, 0]
    top_ys = top_pixel_coordinate_array[::, 1]
    top_zs = top_pixel_coordinate_array[::, 2] 
    bot_xs = bottom_pixel_coordinate_array[::, 0]
    bot_ys = bottom_pixel_coordinate_array[::, 1]
    bot_zs = bottom_pixel_coordinate_array[::, 2] 
            
    #Initialize 3D plot
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection = '3d')
            
        #add substrate points
    ax.scatter(top_xs, top_ys, top_zs, s=.2, color = 'teal', alpha=0.5)
    ax.scatter(bot_xs, bot_ys, bot_zs, s=.2, color = 'brown', alpha=0.5)
            
    if plot_first_layer:
        try:
            ax.scatter(first_points[0], first_points[1], first_points[2], s=10, color = 'blue', alpha=0.1)
        except:
            pass
    if plot_last_layer:
        try:
            ax.scatter(last_points[0], last_points[1], last_points[2], s=10, color = 'purple', alpha=0.1)
        except:
            pass

    def init():
        ax.view_init(elev=10., azim=0)
        return [fig]
            
    def animate(i):
        global azim_step
        ax.view_init(elev=10., azim=i*azim_step)
        print(f"Running {i+1} of {gif_frames}")
        return [fig]
            
    # Animate
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                    frames=gif_frames, interval=gif_integrals, blit=True)
    if gif:
        # Save; Takes a while (~5 minutes?)
        save_format = ".mp4"
        format_name = save_format.replace('.','').upper()
                
        root = Tk()
        save_filename = filedialog.asksaveasfilename(filetypes=[(format_name, save_format)])
        root.destroy()
                
        anim.save(save_filename+save_format, fps=30, extra_args=['-vcodec', 'libx264'])
        #anim.save(save_filename, fps=30)
    else:
        plt.title(voxel_name)
        if save_plot:
            plt.savefig(save_name, dpi=600)
        plt.show()


def voxel_stl_from_gcode(filename_lists = None, 
                        save_stl = True, show_stl = False,
                        flip_normals= True, n_voxel_points = 5,
                        num_strand_exterior_points = 5,
                        minimum_point_dist = 0.2,
                        voxel_size= 0.01, size_multiplier = 1):
    
    """
    Created:   2025-03-26
    Modified:  2025-06-24
    Version:   0.4.0 (beta)

    @author: Aaron Pital (Los Alamos National Lab)

    Description: Taking a 'config.json' file and at least one 'layer.pgm' file with gcode coordinates, take 'n_voxel_points' number of random points and generate
        an STL.
        Wraps up prior versions into a function. Retains cell markers (#%%) for use as a Notebook or in an IDE with appropriate extensions (Spyder natively).

    INPUT:
      (optional)
        filename_lists          2D list with each i,[] list containing a config file and at least one .pgm layer file. If no list is passed, open a dialog to select files.
        save_stl                Save the STL file to the 'filename_lists' directory
        show_stl                Attempt to visualize STL (not currently working)
        flip_normals            Flip normals to point 'out' of strand rather than 'in'
        n_voxel_points          Number of random points to select from parts
        minimum_point_dist      Minimum distance between gcode points allowed before new points are interpolated
        voxel_size              Size (in mm) of each voxel's side-length. Voxels should be 'perfect' cubes, but may be off a bit from rounding errors.
      
      (deprecated; kept for potential future scope)
        num_strand_exterior_points  number of surface points to generate at exterior of of each gcode coordinate 
        size_multiplier         Multiplier to 'blow up' strand size; legacy of prior iterations that's kept for future scope

    OUTPU:
        voxel_dict              Updated with new coordinate measures
    """

    #Hard-coded contraints
    voxel_max_divisor = 10  #voxel can't be more that this size smaller than 'minimum_point_dist'
    voxel_min_divisor = 2   #voxel can't be less that this size smaller than 'minimum_point_dist'
      #random-point voxel parameters
      # Note: different than the voxel array that's returned as an STL; This is the 'camera' voxel that pulls points to feed into the STL voxel
    n_pixels = 100          #number of pixels in the voxel 'camera'; irrelevant for this function, but kept as an option for future scope
    radius = 3              #radius of 'camera' to pull points from (in mm); actually pulls a square region defined by the radius with 2*radius per side

    ##%% (1.2) Select file(s) to make STLs of
    if filename_lists == None:
        continue_check = True
        filename_lists = []
        filenames_idx = 0

        while continue_check:
            # Get the names by user selection
            root = Tk()
            filenames = filedialog.askopenfilenames(title= 'Select images to stitch')
            root.destroy()
    
            # Sort to get images in left-to-right, top-to-bottom order by index (hopefully)
            filename_lists.append(filenames)
            print(f"____Index {filenames_idx}_________")
            for idx, filename in enumerate(filenames):
                print(f"{idx} \t {filename}")
    
            # Ask if user want to add more
            user_report = input("Stitch another set of images (y/n or any key+ENTER to quit)").lower()
            if 'y' in user_report:
                continue_check = True
            else:
                continue_check = False
        
            filenames_idx += 1


    ##%% (2.a) Process file(s)
    for filenames_idx, filenames in enumerate(filename_lists):
    
        # Open and generate a LayerOver part
        try:
            part = Gcode()
            part.get_gcode(files = filenames)
        
                #get names of layers
            layer_list = list(part.layers.keys())
        
                #auto-generate part name guess and assign directory to save to
            part_dir = os.path.dirname(part.config_path)
            part_name_guess = os.path.basename(part_dir)
        
                #check if there's a single strand diameter
            if part.strand_diameter != 0:
                strand_radius = part.strand_diameter/2*size_multiplier
                homogenous_strand_size = True
            else:
                homogenous_strand_size = False
                radius_list = [size/2*size_multiplier for size in part.layer_strand_diameters]


        except Exception as e:
            print(traceback.format_exception(*sys.exc_info()))
            print()
            print("Random failure")
        
        
    ##%%  (2.b) Get random points to evaluate STLs on
    point_dicts = get_random_point_summaries(part, n_retrieval_points = n_voxel_points, \
                                        excluded_point_indcs = [], show_iterations = True)
        
    ##%%  (3.a)
    #Get STL from each point in point_dict
      #grab appropriate keys to index data from dicts
    point_keys = list(point_dicts.keys())
    layer_keys = list(part.layers.keys())

    for point_key in point_keys:
    
        try:
            for layer_idx, layer_key in enumerate(layer_keys):
                
                #generate voxel name
                voxel_name = str(part_name_guess + "_voxel-" + str(layer_idx))
                voxel_point_name = str(part_name_guess + "_VoxelPoints" + str(layer_idx))
                save_name = os.path.join(part_dir, voxel_point_name)
                
                #pull voxel features
                center = point_dicts[layer_idx]['approximate_center']
                normal_vector = point_dicts[layer_idx]['unit_normal']
                #TODO: Radius is currently hard-coded by user or default input; could be made to be more flexible here
                
                #get voxel bounds and add to dictionary 
                voxel_dict = get_voxel(center, normal_vector, radius= radius, voxel_depth= 2*radius, n_pixels=n_pixels)
                voxel_dict = get_layer_points(part, voxel_dict, boundary_buffer_multiplier = 0.2)

                #update voxel_dict
                voxel_dict.update({'layer_names': layer_list})
                voxel_dict.update({'voxel_name': voxel_name})
                voxel_dict.update({'save_filepath': save_name})
                
                #Visualize voxel
                voxel_points_visualize(voxel_dict, 
                           gif = False,
                           gif_frames = 4,
                           gif_integrals = 50,
                           azim_total_degrees = 360,
                           azim_step = 0,
                           save_plot = True)
                    
                #Flatten array of coordinates to make it homogenous
                #   (otherwise it's a list of lists, and each sub-list is a different size; this makes numpy sad)
                voxel_dict.update( {'flat_layer_coordinates':{} } )
                for layer_key in list(part.layers.keys()):
                    #flatten coordinate array; kept previously as inhomogenous list to make it easier to parse consecutive strand coordinates
                       #NOTE: each 'this_list' is a strand of the print visible within the 'camera voxel' 
                    flat_coordinates = []
                    for this_list in voxel_dict['layer_coordinates'][layer_key]['coordinates']:
                        for point in this_list:
                            flat_coordinates.append(point)
                
                    #re-wrap lists as numpy arrays to make the following stuff faster and the syntax clearer (mostly the faster part)
                    flat_coordinates = np.array(flat_coordinates)
                    flat_indices = np.array(voxel_dict['layer_coordinates'][layer_key]['indices'])
                    
                    #write the flattened results to the returned dict
                    voxel_dict['flat_layer_coordinates'].update ( {layer_key: {
                                            'coordinates':flat_coordinates, \
                                            'indices':flat_indices} })
                
                #Define variables for this layer
                this_save_name = f"{part_name_guess}_{layer_key.replace('.txt','')}"
                layer_nozzle_size = part.layer_strand_diameters[layer_idx]
                strand_radius = layer_nozzle_size/2

                    #Pull the flat coordinates and their indices (made just above)
                these_coordinates = voxel_dict['flat_layer_coordinates'][layer_key]['coordinates']
                these_indices = voxel_dict['flat_layer_coordinates'][layer_key]['indices']

                    # Get max/min boundaries for setting up smart grid coordinates
                max_x = np.max(these_coordinates[::,0])
                min_x = np.min(these_coordinates[::,0])
                max_y = np.max(these_coordinates[::,1])
                min_y = np.min(these_coordinates[::,1])
                max_z = np.max(these_coordinates[::,2])
                min_z = np.min(these_coordinates[::,2])
                x_dim = abs(max_x-min_x)
                y_dim = abs(max_y-min_y)
                z_dim = abs(max_z-min_z)
                # Make sure the dimensions allow for full strand size
                if x_dim < (strand_radius*2)*1.5:
                    x_dim = (strand_radius*2)*1.5
                if y_dim < (strand_radius*2)*1.5:
                    y_dim = (strand_radius*2)*1.5
                if z_dim < (strand_radius*2)*1.5:
                    z_dim = (strand_radius*2)*1.5
                
                    # adjust parameters for new size constraints
                if (layer_nozzle_size/minimum_point_dist < 2):
                    minimum_point_dist = layer_nozzle_size/2
                    print(f"Adjusting 'minimum_point_dist' to: {round(minimum_point_dist, 3)}")
                if (minimum_point_dist/voxel_size > voxel_max_divisor) or (minimum_point_dist/voxel_size < voxel_min_divisor):
                    #Make sure the voxel size is 'well-matched' to strand diameter
                    voxel_size = round(minimum_point_dist/voxel_min_divisor, 5)
                voxel_x_grid_size = round(x_dim/voxel_size) + 10
                voxel_y_grid_size = round(y_dim/voxel_size) + 10
                voxel_z_grid_size = round(z_dim/voxel_size) + 10

                # Convert a point cloud to a voxel grid; Create an empty voxel grid
                    #create the global voxel array
                voxels = np.zeros((voxel_x_grid_size, voxel_y_grid_size, voxel_z_grid_size))
                
                    #create arrays of actual locations to index point coordinates to voxel array indices
                half_voxel_size = voxel_size/2
                voxel_x_locations = np.linspace(min_x-(5*half_voxel_size), max_x+(5*half_voxel_size), voxel_x_grid_size)
                voxel_y_locations = np.linspace(min_y-(5*half_voxel_size), max_y+(5*half_voxel_size), voxel_y_grid_size)
                voxel_z_locations = np.linspace(min_z-(5*half_voxel_size), max_z+(5*half_voxel_size), voxel_z_grid_size)
                   #generate meshed grid of locations (3, x_dim, y_dim, z_dim)
                mesh_locations = np.array(np.meshgrid(voxel_x_locations,voxel_x_locations,voxel_z_locations, indexing='ij'))
                   #populate voxel_locations array
                   #TODO: this is inelegant, but I'm having trouble getting my head around an alternative
                voxel_locations = np.zeros((voxel_x_grid_size, voxel_y_grid_size, voxel_z_grid_size, 3))
                voxel_locations[::, ::, ::, 0] = mesh_locations[0, ::, ::, ::]
                voxel_locations[::, ::, ::, 1] = mesh_locations[1, ::, ::, ::]
                voxel_locations[::, ::, ::, 2] = mesh_locations[2, ::, ::, ::]
            
                #Initialize variables for looping through gcode points
                points = []
                normals = []
                    # initialize 'last_position' to first coordinate and get first radial (strand exterior) points
                last_position = these_coordinates[0, ::]
                last_index = these_indices[0]
                #Get coordinates and normals
                pbar = tqdm(total= part.layers[layer_key]['coordinates'].shape[0])  #initialize progress bar
                for current_index, current_position in zip(these_indices[1::], these_coordinates[1::,::]):
                    pbar.update(1)  #update progress bar
                    points.append(current_position)
                    
                    #Check if this point is in the same strand as the last point. If not, move on.
                    if (current_index-last_index)>1:
                        continuous_point = False

                    else:
                        continuous_point = True
                        #Get this point and the 'last point', calculate max distance bewteen strand surface points
                            # max distance is bottom of one strand at 'p1' to top of strand at 'p2'; calc. that max distance
                        p1 = np.array([last_position[0], last_position[1], 0])
                        p2 = np.array([current_position[0], current_position[1], layer_nozzle_size])
                        distance = math.sqrt( (p1[0]-p2[0])**2 + \
                                                (p1[1]-p2[1])**2 + \
                                                (p1[2]-p2[2])**2)

                    #Run point-by-point, check if next point is <= 'minimum_point_distance', and interpolate new points if it's not
                    if continuous_point and (distance > minimum_point_dist):
                        num_points = int(distance / minimum_point_dist) + 1
                        if num_points <2:
                            num_points = 2
                        interpolated_points = [\
                                (last_position[0] + (current_position[0] - last_position[0]) * t/num_points, \
                                    last_position[1] + (current_position[1] - last_position[1]) * t/num_points, \
                                    last_position[2] + (current_position[2] - last_position[2]) * t/num_points) for t in range(1, num_points+1)]
                    
                        for current_position in interpolated_points:
                            #Calculate direction vector
                            dx = current_position[0] - last_position[0]
                            dy = current_position[1] - last_position[1]
                            dz = current_position[2] - last_position[2]
                            if flip_normals:
                                #rotate -90 degrees to get outward normal (clockwise)
                                dx, dy = dy, -dx
                            #normalize the normal vector
                            length = math.sqrt(dx**2 + dy**2 + dz**2)
                            dx = dx/length
                            dy = dy/length
                            dz = dz/length
                            point_normal_vector = [dx, dy, dz]
                            normals.append(point_normal_vector)

                            #Get voxel to line-segment (strand) distance matrix
                                #focus in on only voxel region this segment occupies
                            vox_x_min = min(current_position[0], last_position[0])
                            vox_x_max = max(current_position[0], last_position[0])
                            vox_y_min = min(current_position[1], last_position[1])
                            vox_y_max = max(current_position[1], last_position[1])
                            vox_z_min = min(current_position[2], last_position[2])
                            vox_z_max = max(current_position[2], last_position[2])
                            x_min_idx = np.argmin(np.abs(voxel_x_locations - vox_x_min))
                            x_max_idx = np.argmin(np.abs(voxel_x_locations - vox_x_max))
                            y_min_idx = np.argmin(np.abs(voxel_y_locations - vox_y_min))
                            y_max_idx = np.argmin(np.abs(voxel_y_locations - vox_y_max))
                            z_min_idx = np.argmin(np.abs(voxel_z_locations - vox_z_min))
                            z_max_idx = np.argmin(np.abs(voxel_z_locations - vox_z_max))
                                #make sure indices are actually the max or min in absolute value
                            x_min_idx = min(x_min_idx, x_max_idx)
                            x_max_idx = max(x_min_idx, x_max_idx)
                            y_min_idx = min(y_min_idx, y_max_idx)
                            y_max_idx = max(y_min_idx, y_max_idx)
                            z_min_idx = min(z_min_idx, z_max_idx)
                            z_max_idx = max(z_min_idx, z_max_idx)
                                #pull sub-region of interest
                            sub_voxel = voxels[x_min_idx:x_max_idx,
                                                y_min_idx:y_max_idx,
                                                z_min_idx:z_max_idx]
                            sub_voxel_locations = voxel_locations[x_min_idx:x_max_idx,
                                                y_min_idx:y_max_idx,
                                                z_min_idx:z_max_idx,
                                                ::]
                            sub_voxel_array = sub_voxel_locations.reshape(3,-1).T  #make the m,n,o,3 array into an mxnxo,3 array for speed
                                #ironically, get rid of that speed by adding a FOR loop
                            #TODO: vectorize and speed this up by a lot
                            distance_array = []
                            for idx in range(sub_voxel_array.shape[0]):
                                distance_array.append(Points.point_line_distance(sub_voxel_array[idx], last_position, current_position))
                            distance_array = np.array(distance_array)
                                #shape back into orginal size 
                            distance_voxel = distance_array.reshape(x_max_idx-x_min_idx,
                                                                    y_max_idx-y_min_idx,
                                                                    z_max_idx-z_min_idx)
                                #mark all of the global 'voxels' voxels that are "within" the strand as calculated above
                            sub_voxel[distance_voxel <= strand_radius] = 1
                                
                            last_position = current_position
                    
                    #If this point is continous and the distance to the last point doesn't require interpolation, just do the math
                    elif continuous_point:
                        #Calculate direction vector
                        dx = current_position[0] - last_position[0]
                        dy = current_position[1] - last_position[1]
                        dz = current_position[2] - last_position[2]
                        if flip_normals:
                            #rotate -90 degrees to get outward normal (clockwise)
                            dx, dy = dy, -dx
                        #normalize the normal vector
                        length = math.sqrt(dx**2 + dy**2 + dz**2)
                        dx = dx/length
                        dy = dy/length
                        dz = dz/length
                        point_normal_vector = [dx, dy, dz]
                        normals.append(point_normal_vector)

                        #Get voxel to line-segment (strand) distance matrix
                            #focus in on only voxel region this segment occupies
                        vox_x_min = min(current_position[0], last_position[0])
                        vox_x_max = max(current_position[0], last_position[0])
                        vox_y_min = min(current_position[1], last_position[1])
                        vox_y_max = max(current_position[1], last_position[1])
                        vox_z_min = min(current_position[2], last_position[2])
                        vox_z_max = max(current_position[2], last_position[2])
                        x_min_idx = np.argmin(np.abs(voxel_x_locations - vox_x_min))
                        x_max_idx = np.argmin(np.abs(voxel_x_locations - vox_x_max))
                        y_min_idx = np.argmin(np.abs(voxel_y_locations - vox_y_min))
                        y_max_idx = np.argmin(np.abs(voxel_y_locations - vox_y_max))
                        z_min_idx = np.argmin(np.abs(voxel_z_locations - vox_z_min))
                        z_max_idx = np.argmin(np.abs(voxel_z_locations - vox_z_max))
                            #make sure indices are actually the max or min in absolute value
                        x_min_idx = min(x_min_idx, x_max_idx)
                        x_max_idx = max(x_min_idx, x_max_idx)
                        y_min_idx = min(y_min_idx, y_max_idx)
                        y_max_idx = max(y_min_idx, y_max_idx)
                        z_min_idx = min(z_min_idx, z_max_idx)
                        z_max_idx = max(z_min_idx, z_max_idx)
                            #pull sub-region of interest
                        sub_voxel = voxels[x_min_idx:x_max_idx,
                                            y_min_idx:y_max_idx,
                                            z_min_idx:z_max_idx]
                        sub_voxel_locations = voxel_locations[x_min_idx:x_max_idx,
                                            y_min_idx:y_max_idx,
                                            z_min_idx:z_max_idx, 
                                            ::]
                        sub_voxel_array = sub_voxel_locations.reshape(3,-1).T  #make the m,n,o,3 array into an mxnxo,3 array for speed
                        
                        #ironically, get rid of that speed by adding a FOR loop
                          #TODO: vectorize and speed this up by a lot
                        distance_array = []
                        for idx in range(sub_voxel_array.shape[0]):
                            distance_array.append(Points.point_line_distance(sub_voxel_array[idx], last_position, current_position))
                            #shape back into orginal size 
                        distance_voxel = distance_array.reshape(x_max_idx-x_min_idx,
                                                                y_max_idx-y_min_idx,
                                                                z_max_idx-z_min_idx)
                            #mark all of the global 'voxels' voxels that are "within" the strand as calculated above
                        sub_voxel[distance_voxel <= strand_radius] = 1

                        last_position = current_position

                    else:
                        last_position = current_position
                        normals.append([0,0,0])
                    
                    
                pbar.close()  #close progress bar
            
            #Save as .xyz file
            xyz_save_name = this_save_name + '.xyz'
            xyz_save_filepath = os.path.join(part_dir, xyz_save_name)
            with open(xyz_save_filepath, 'w') as file:
                for point, normal in zip(points, normals):
                    file.write(f"{point[0]} {point[1]} {point[2]} {normal[0]} {normal[1]} {normal[2]}\n")
              
        
            #Convert a voxel grid to an STL file using the Marching Cubes algorithm, ensuring the output matches the original scale.
                # calculate the scale factors based on the point cloud dimensions and the voxel grid
            min_bound = np.min(points, axis=0)
            max_bound = np.max(points, axis=0)
            scales = (max_bound - min_bound) / np.array(voxels.shape)
    
            verts, faces, _, _ = marching_cubes(voxels)
        
            #     # scale vertices back to the original point cloud dimensions
            # verts = verts * scales + min_bound
            stl_mesh = Mesh(np.zeros(faces.shape[0], dtype= Mesh.dtype))
            if flip_normals:
                # Reverse the order of vertices for each face to flip normals
                faces = faces[:, ::-1]
            for i, f in enumerate(faces):
                for j in range(3):
                    stl_mesh.vectors[i][j] = verts[f[j], :]
        
            #Save/show the results
                # get save name and filepath
            stl_save_name = this_save_name + '.stl'
            stl_save_filepath = os.path.join(part_dir, stl_save_name)
        
            if save_stl:
                stl_mesh.save(stl_save_filepath)
        
            if show_stl:
                #Open3D approach
                # vis = o3d.visualization.Visualizer()
                # vis.create_window()
                # vis.add_geometry(stl_mesh)
                # vis.run()
                # vis.destroy_window()
            
                #   #attempt to load stl to ensure saving worked
                # opened_mesh = o3d.io.read_triangle_model(stl_save_filepath)
                # visualize(opened_mesh)
            
                #Matplotlib approach
                # Create a new plot
                figure = plt.figure()
                axes = figure.add_subplot(projection='3d')
            
                # Load the STL files and add the vectors to the plot
                opened_mesh = Mesh.from_file(stl_save_filepath)
                poly_collection = mplot3d.art3d.Poly3DCollection(opened_mesh.vectors)
                poly_collection.set_color((0.7,0.7,0.7))  # play with color
                axes.add_collection3d(poly_collection)
            
                # Show the plot to the screen
                plt.show()
    
        except Exception as e:
            print()
            print('#'*50)
            print(f"Exception encountered: {e}")
            print(traceback.format_exc())
            print('#'*50)
            print()
        
        # |  To-be-deleted...
        # V  Duplicated from code import? Keeping it in case I'm wrong and because I'm too lazy to figure it out right now

        # #Save as .xyz file
        # this_save_name = f"{part_name_guess}_points-normals"
        # xyz_save_name = this_save_name + '.xyz'
        # xyz_save_filepath = os.path.join(part_dir, xyz_save_name)
        # with open(xyz_save_filepath, 'w') as file:
        #     for point, normal in zip(points, normals):
        #         file.write(f"{point[0]} {point[1]} {point[2]} {normal[0]} {normal[1]} {normal[2]}\n")
           
        # # Convert a point cloud to a voxel grid; Create an empty voxel grid
        # voxels = np.zeros((voxel_x_grid_size, voxel_y_grid_size, voxel_z_grid_size))
        # min_bound = np.min(points, axis=0) - voxel_size
        # max_bound = np.max(points, axis=0) + voxel_size
        # scales = (max_bound - min_bound) / voxel_size
        # indices = ((points - min_bound) / scales).astype(int)
        # for index in indices:
        #     voxels[index[0], index[1], index[2]] = 1
    
        # #Convert a voxel grid to an STL file using the Marching Cubes algorithm, ensuring the output matches the original scale.
        #     # calculate the scale factors based on the point cloud dimensions and the voxel grid
        # min_bound = np.min(points, axis=0)
        # max_bound = np.max(points, axis=0)
        # scales = (max_bound - min_bound) / np.array(voxels.shape)

        # verts, faces, _, _ = marching_cubes(voxels)
    
        #     # scale vertices back to the original point cloud dimensions
        # verts = verts * scales + min_bound
        # stl_mesh = Mesh(np.zeros(faces.shape[0], dtype= Mesh.dtype))
        # if flip_normals:
        #     # Reverse the order of vertices for each face to flip normals
        #     faces = faces[:, ::-1]
        # for i, f in enumerate(faces):
        #     for j in range(3):
        #         stl_mesh.vectors[i][j] = verts[f[j], :]
    
        # #Save/show the results
        #     # get save name and filepath
        # stl_save_name = this_save_name + '.stl'
        # stl_save_filepath = os.path.join(part_dir, stl_save_name)
    
        # if save_stl:
        #     stl_mesh.save(stl_save_filepath)
    
        # if show_stl:
        #     #Open3D approach
        #     # vis = o3d.visualization.Visualizer()
        #     # vis.create_window()
        #     # vis.add_geometry(stl_mesh)
        #     # vis.run()
        #     # vis.destroy_window()
        
        #     #   #attempt to load stl to ensure saving worked
        #     # opened_mesh = o3d.io.read_triangle_model(stl_save_filepath)
        #     # visualize(opened_mesh)
        
        #     #Matplotlib approach
        #     # Create a new plot
        #     figure = plt.figure()
        #     axes = figure.add_subplot(projection='3d')
        
        #     # Load the STL files and add the vectors to the plot
        #     opened_mesh = Mesh.from_file(stl_save_filepath)
        #     poly_collection = mplot3d.art3d.Poly3DCollection(opened_mesh.vectors)
        #     poly_collection.set_color((0.7,0.7,0.7))  # play with color
        #     axes.add_collection3d(poly_collection)
        
        #     # Show the plot to the screen
        #     plt.show()

    return voxel_dict