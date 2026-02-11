# -*- coding: utf-8 -*-
"""
? 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2024-03-18
Modified:  2026-02-11
Version:   0.6.1

@author: Aaron Pital (Los Alamos National Lab)

Description: Class wrapper for handling 3D data; parametric conversion, 3D transformations, 
    point interpolation and convex hull generation, and point/object visualization.
    
    Updates:
        2025-03-27:
            - Added random initialization support to 'generate_radial_points'
        2025-04-16
            - Fixed # of points returned in 'generate_radial_points'
        2025-07-15
            - Refactored STL
        2026-02

TODO:
    get_radial_neighbor_array_data
        -Add functionality for interior and mixed interior/edge points

"""


  #System and built-ins
import math
import random
from tkinter import Tk, filedialog

  #Visualizaiton
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import animation

  #Data Handling
import numpy as np
from skimage.draw import line

  #Scientifiic algorithm packages
from scipy.interpolate import interp1d
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

    
#take G-Code coordinates parsed from file and return a layer dictionary with XYZ-cartesian coordinates
def AC_to_XYZ(xyzacf_array, substrate_dict, offset_adjust=False):
    ''' v1.0   created:2024-03-22   modified:2024-03-24
    Take a n,6 sized array (6) of xyz-AzimuthalAngle(a)-RotationAngle(c) array and convert each row to a set of XYZ coordinates in Cartesian space.
    As of v1.0, only coordinates are logged. Future versions could use feedrate and some fancier equations to produce.
    INPUT:  'xyzacf_array'- np.array(n, [x,y,z,a,c,f,mode])
            'substrate_dict'- <gcode_object>.substrate_dict
    
    OUTPUT: 'layer_dict'
    '''
        
    # Initialize variables
    layer_dict = {
        'coordinates':     np.array([]),
        'flat_part':       False,
        'global_feedrate': 0
    }
    Xs = []
    Ys = []
    Zs = []
    devl_testing = []
        
    # Parse coordinates; lowercase parameters are trial, polar coordinates until conversion
    for point in xyzacf_array:
        x, y, z, a, c, f, mode = point
            
        # If part is flat, first line should be blank except for index 5 and mode='F'
        if mode[0]=='F':
            layer_dict['global_feedrate'] = f
            layer_dict['flat_part'] = True
        
        # Parse simple 'G1' coordinate moves; 
        if mode=='G1':
            # Parse coordinates directly if part is flat
            if layer_dict['flat_part']:
                Xs.append(x)
                Ys.append(y)
                Zs.append(z)
                    
            # If part is not flat and requires offset adjustments, parse xyzac coordinates to XYZ cartesian coordinates and apply offsets
            elif offset_adjust:
                #TODO: create framework to parse and apply offsets
                pass
        
            # If part is not flat and no offsets required, directly parse xyzac coordinates to XYZ cartesian coordinates
            else:
                #convert the A-C-Z coordinates to XYZ
                c_radians = math.pi*c/180
                a_radians = math.pi*a/180
                #get level-curve radius at angle
                    #find index of surface profile closest to 'a'
                min_idx = np.nanargmin(np.abs(substrate_dict['surface_profile'][:,0]- a))  #idx 0-angle, 1-lat, 2-vert
                this_lat = substrate_dict['surface_profile'][min_idx,1]  # radius from polar centerline; 0@ 0, full radius@ 90
                this_height = substrate_dict['surface_profile'][min_idx,2]
                    
                #trig X-Y; used to calculate Z but must be adjusted by x and y to get cartesian X and Y coordinates
                acX = math.cos(c_radians) * this_lat  #'z' is radius value
                acY = math.sin(c_radians) * this_lat
                    #adjust to final X and Y
                X = acX + x
                Y = acY + y
                    
                #trig out the Z; right triangle of angle 'a' and hypoteneuse 'z'
                    # Attempts that failed because they don't follow level-curve @ 'vert' of mandrel
                    #adjacent_leg = math.sqrt( (acX**2 + acY**2) )
                    #adjacent_leg = math.sqrt( (X**2 + Y**2) )
                    #Z = adjacent_leg * tan_a   #failed because tangent is a shitty function
                      
                    #this should be sin()... not entirely sure why this works tbh
                    #TODO: figure out why this works and why future stuff will probably break here
                Z = (z) * math.cos(a_radians)
        
                #write the final cartesian coordinates
                Xs.append(X)
                Ys.append(Y)
                Zs.append(Z)
        
    XYZ_array = np.stack((Xs,Ys,Zs), axis=1)
    layer_dict['coordinates'] = XYZ_array
        
    return layer_dict
    

#Simple unit-normal function for getting points on a circle
def circle_circumference_xy(normal_vector, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    normal_vector = np.asarray(normal_vector)
    normal_vector = normal_vector / math.sqrt(np.dot(normal_vector, normal_vector))
    a = math.cos(theta / 2.0)
    b, c, d = -normal_vector * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    
#Generate points on circle circ. given line, reference point, and radius
def generate_radial_points(vector, reference_point, radius, n_circ_points = 6, randomize_radial_start = False):
    ''' v0.2.3 created:2024-05-12  modified:2025-04-17 
        
    Simple function to get 'n_circ_points' number of circumferential points at 'radius' distance from 'reference_point' along
        the line 'normal_vector'.
          
    Updates:
        - 2025-03-37: added 'randomize_radial_start' boolean; adds a random start radial angle
        - 2025-04-16: added 'angle_step_size' to fix number of points returned
        
    '''
        
    #Get angle step to obtain given # of points
    angle_step_size = 360//n_circ_points
        
    #Make sure we're numpying 
    if type(vector) == 'list':
        vector = np.array(vector)
            
    #Get unit vector and ascribe to radius
    x = np.array([1,0,0]).astype(np.float64) # take a random vector of magnitude 1
    x -= np.multiply(x.dot(vector), vector) / np.linalg.norm(vector)**2  # make it orthogonal to n
    x /= np.linalg.norm(x)  # normalize
        #find first point on circle (x1). 
        #currently it has magnitude of 1, so we multiply it by the r
    x1 = reference_point + (x*radius)
    
        #generate vector from centre ('reference_point') to first circle point
    center_circumference_vector = x1 - reference_point
    
    #Run through points on a circle and spit out coordinates
        # apply random '0' angle start if called for
    if randomize_radial_start:
        zero_angle = random.randint(0, 360)
    else:
        zero_angle = 0
        # generate the actual points
    cirlce_circumference_points = []
    for theta in range(0, 360, angle_step_size):
        theta = zero_angle + theta
        while theta > 360:
            if theta > 360:
                theta = theta-360
        circle_i = np.dot(circle_circumference_xy(vector, np.deg2rad(theta)), center_circumference_vector)
        point = circle_i+reference_point
        cirlce_circumference_points.append(point)
    
    return cirlce_circumference_points
    
    
# Takes in spherical coordinates and makes a convex hull
def convexHull(xcoords, ycoords, zcoords, \
                center = [0,0,0]):
    
    ''' v1.1  created:2024-03-19  modified:2024-03-21
    (Modified from https://stackoverflow.com/questions/77564155/how-can-one-plot-a-3d-surface-in-matplotlib-by-points-coordinates)
    Take points, project to a sphere (sphere-ify), take convex hull, and shrink back down to original scale (un-sphere-ify)    
    '''
    # Adjust for 2D flat plate 
    flat = False
    unique_zs = np.unique(np.array(zcoords))
    if (sum(np.unique(np.array(zcoords)))<1.02 and sum(np.unique(np.array(zcoords)))>0.98) or \
        (sum(np.unique(np.array(zcoords)))<0.02 and sum(np.unique(np.array(zcoords)))>-0.02) :
        zcoords = np.ones(len(zcoords))
        flat = True
        
    # Combine points into a single array
    original_pts = np.stack((xcoords,ycoords,zcoords), axis=1)
    
    # Generate point data
        #calculate the location of each point when it is expanded out to the sphere
    kdtree = KDTree(original_pts) # tree of nearest points
        #'d' is an array of distances, 'sphereInd' is array of indices
    d, sphereIndcs = kdtree.query(center, original_pts.shape[0])
    spherePts = np.zeros(original_pts.shape, dtype=float)
        #sphere-ify; expand to sphere
    radius = np.amax(d)
    for p in range(original_pts.shape[0]):
        spherePts[p] = original_pts[sphereIndcs[p]] *radius /d[p]
    
    # Generate convex hull and simplices
    if not flat:
        # For mandrels that are spherical-ish curves
        hull = ConvexHull(spherePts)
        triangInds = hull.simplices # returns the list of indices for each triangle
            #clean big points
        triangInds = removeBigTriangs(original_pts[sphereIndcs], triangInds)
        return original_pts, sphereIndcs, triangInds
            
    else:
        # For flat plates
        twoD_pts = np.stack((xcoords,ycoords), axis=1)
        hull = Delaunay(twoD_pts)
        triangInds = hull.simplices # returns the list of indices for each triangle
            #clean big points
        #triangInds = point_clod.removeBigTriangs(original_pts[sphereIndcs], triangInds)
        return original_pts, sphereIndcs, triangInds
            

def check_point_line_distance(point, line_vector, reference_points = np.array([0,0,0]), reference_delta=0):
    ''' v0.2.0  created:2024-05-12  modified:2024-05-20 
        
    INPUT:   Point, line vector, and location bounds on line to define distance region
                'reference_points'- np ndarray of either 1 or two points. If one, supply a 'reference_delta' to ge a line.
                    If two points supplied, make a vector between them 
                'reference_delta'- distance above and below reference along 'line' vector to generate reference points
    ACTION:  Calculate distance between line and point and reference_point and point
    OUTPUT:  point_line_distance, point_reference_distance
    '''
    
    #Make sure the reference points
        
    #If there's a 'reference_delta', assume there's a reference point too and do everything wrt those
        #assume vector is normal
    if reference_delta != 0 and len(reference_points.shape) == 1:
        #Generate bounding points for centerline
        rp_lower = reference_points - (reference_delta * line_vector)
        rp_upper = reference_points + (reference_delta * line_vector)
        
        #normalized tangent vector
        d = np.divide(rp_lower - rp_upper, np.linalg.norm(rp_lower - rp_upper))
        
        #signed parallel distance components
        s = np.dot(rp_upper - point, d)
        t = np.dot(point - rp_lower, d)
        
        #clamped parallel distance
        h = np.maximum.reduce([s, t, 0])
            
        #perpendicular distance component
        c = np.cross(point - rp_upper, d)
    
    #Otherwise (most general case), just consider point -> line distance without assuming anything
    else:
        p1,p2 = reference_points
        
        #normalized tangent vector
        d = np.divide(p2-p1, np.linalg.norm(p2-p1))
        
        #signed parallel distance components
        s = np.dot(p2 - point, d)
        t = np.dot(point - p1, d)
        
        #clamped parallel distance
        h = np.maximum.reduce([s, t, 0])
            
        #perpendicular distance component
        c = np.cross(point - p2, d)
        
    return np.hypot(h, np.linalg.norm(c))


def line_line_distance(a0,a1,b0,b1,\
                                clampAll=False,clampA0=False,clampA1=False,clampB0=False,clampB1=False):
    ''' v0.1.0  created:2024-05-20  modified:2024-05-20 
    taken verbatim from https://stackoverflow.com/questions/2824478/shortest-distance-between-two-line-segments
        
    INPUT:   Two 'a' and 'b' points; point a/b arbitrary, as is 0/1
                'Clamp' options constrain distances within those segment bound(s); otherwise, the shortest line distance is calculated.
    ACTION:  lorem
    OUTPUT:  lorem
    '''
    
    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0=True
        clampA1=True
        clampB0=True
        clampB1=True
    
    
    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0
    magA = np.linalg.norm(A)
    magB = np.linalg.norm(B)
        
    _A = A / magA
    _B = B / magB
        
    cross = np.cross(_A, _B);
    denom = np.linalg.norm(cross)**2
        
        
    # If lines are parallel (denom=0) test if lines overlap.
    # If they don't overlap then there is a closest point solution.
    # If they do overlap, there are infinite closest positions, but there is a closest distance
    if not denom:
        d0 = np.dot(_A,(b0-a0))
            
        # Overlap only possible with clamping
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A,(b1-a0))
                
            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 and clampB1:
                    if np.absolute(d0) < np.absolute(d1):
                        return a0,b0,np.linalg.norm(a0-b0)
                    return a0,b1,np.linalg.norm(a0-b1)
                    
                    
            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return a1,b0,np.linalg.norm(a1-b0)
                    return a1,b1,np.linalg.norm(a1-b1)
                    
                    
        # Segments overlap, return distance between parallel segments
        return None,None,np.linalg.norm(((d0*_A)+a0)-b0)
            
        
        
    # Lines criss-cross: Calculate the projected closest points
    t = (b0 - a0);
    detA = np.linalg.det([t, _B, cross])
    detB = np.linalg.det([t, _A, cross])
    
    t0 = detA/denom;
    t1 = detB/denom;
    
    pA = a0 + (_A * t0) # Projected closest point on segment A
    pB = b0 + (_B * t1) # Projected closest point on segment B
    
    
    # Clamp projections
    if clampA0 or clampA1 or clampB0 or clampB1:
        if clampA0 and t0 < 0:
            pA = a0
        elif clampA1 and t0 > magA:
            pA = a1
            
        if clampB0 and t1 < 0:
            pB = b0
        elif clampB1 and t1 > magB:
            pB = b1
                
        # Clamp projection A
        if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
            dot = np.dot(_B,(pA-b0))
            if clampB0 and dot < 0:
                dot = 0
            elif clampB1 and dot > magB:
                dot = magB
            pB = b0 + (_B * dot)
        
        # Clamp projection B
        if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
            dot = np.dot(_A,(pB-a0))
            if clampA0 and dot < 0:
                dot = 0
            elif clampA1 and dot > magA:
                dot = magA
            pA = a0 + (_A * dot)
    
        
    return pA,pB,np.linalg.norm(pA-pB)
        
#
def fibonnaci_points(lats, verts, points=1000):
    if points <= 1:
        points = len(verts)
    else:
        if points > len(verts):
            lats, verts = curve_interpolate(lats, verts, points)
                
    # Make some golden ratio
    phi = math.pi*(math.sqrt(5)-1)

    # For each 'y' provided in 'verts' 
    xs=[]
    ys=[]
    zs=[]
    radii = []
    for i in range(points):
        theta = phi*i
        y = verts[i]
        this_radius = lats[i]
        x = math.cos(theta)*this_radius
        z = math.sin(theta)*this_radius
        
        xs.append(x)
        ys.append(y)
        zs.append(z)
        radii.append(math.sqrt(x**2+y**2+z**2))
    return xs, ys, zs, radii


def get_3point_normal(point_coords, show_plot = False):

    ''' v1.2.0  created:2024-05-01  modified:2024-05-12
        Take a set of three x,y,z points and return a dictionary with those points, their normal vector, and other stuff.
    '''
        
    return_dict = {}
    return_dict.update({'original_points': point_coords})
        
    # Get normal vector and point lists
    xs = []
    ys = []
    zs = []
    for point in point_coords:
        xs.append(point[0])
        ys.append(point[1])
        zs.append(point[2])

    #Extract points to 'pX's, define cartesian coordinates
    p0, p1, p2 = point_coords
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    #Generate direction vectors along plane
    ux, uy, uz = u = [x1-x0, y1-y0, z1-z0] #first vector
    vx, vy, vz = v = [x2-x0, y2-y0, z2-z0] #sec vector

    #Take cross product to get normal to plane vectors
    u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx] #cross product
        
    #Get center point, make normal unit vector, add everything to dictionary
    center  = np.array([np.average(xs), np.average(ys), np.average(zs)])   #Roughly middle of 3 points; can be off
    return_dict.update({'center_plane_vector': np.array([x0-center[0], y0-center[1], z0-center[0]])})
    raw_normal = np.array(u_cross_v)
    normal = raw_normal/np.linalg.norm(raw_normal)
    return_dict.update({'approximate_center': center})
    return_dict.update({'raw_normal': raw_normal})
    return_dict.update({'unit_normal': normal})
        
    #Check normal is actually normal; dot product should be ~zero if so
    d = (center-p0).dot(normal)
    return_dict.update({'dot_product': d})

    #Guess at a normal vector size based on triangle-point spaceing for better visualization
    center_point_distances = []
    for point in [p0, p1, p2]:
        xp, yp, zp = point
        this_dist = math.sqrt((center[0]-xp)**2  + (center[1]-yp)**2 + (center[2]-zp)**2)
        center_point_distances.append(this_dist)
    max_distance = max(center_point_distances)/4
    scaled_norm = max_distance * normal
    return_dict.update({'scaled_norm': scaled_norm})

    #Calculate angle between points
    vector_1 = p0-p1
    vector_2 = p0-p2
    this_angle = math.acos( ((vector_1[0]*vector_2[0]+vector_1[1]*vector_2[1]+vector_1[2]*vector_2[2])/  \
            (math.sqrt(vector_1[0]**2+vector_1[1]**2+vector_1[2]**2) * \
                math.sqrt(vector_2[0]**2+vector_2[1]**2+vector_2[2]**2))))
    return_dict.update({'p0_angle': this_angle})
    
    if show_plot:
        # plot the surface
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111, projection = '3d')
        ax.scatter(xs, ys, zs, s=200)
        #ax.scatter(normal[0], normal[1], normal[2], s=200, color = 'red')
        ax.scatter(np.average(xs), np.average(ys), np.average(zs), s=200, color = 'aquamarine')
        ax.quiver(center[0], center[1], center[2], \
                    scaled_norm[0], scaled_norm[1], scaled_norm[2], color="m", \
                    arrow_length_ratio = 0.3, linewidth = 5)
        ax.plot_trisurf(xs, ys, zs, color = 'gray', alpha = 0.5)
        plt.show()

    return return_dict

    
#Helper function for culling invalid polygons
def removeBigTriangs(points, indices, tolerance=10):
    ''' v1.0  created:2024-03-19  modified:2024-03-19
    (Modified from https://stackoverflow.com/questions/77564155/how-can-one-plot-a-3d-surface-in-matplotlib-by-points-coordinates)
    Calculate Euclidean distance and cull points that are outside tolerance range
    '''
    newInds = []
    for idx in indices:
        if ((np.sqrt(np.sum((points[idx[0]]-points[idx[1]])**2, axis=0))<tolerance) and
            (np.sqrt(np.sum((points[idx[0]]-points[idx[2]])**2, axis=0))<tolerance) and
            (np.sqrt(np.sum((points[idx[1]]-points[idx[2]])**2, axis=0))<tolerance)):
            newInds.append(idx)
    return np.array(newInds)
        

# Helper function for visualization and animation
def visualize_3D(xyzarray, sphIndices, triangIndices, \
                    point_color = 'w', map_color = cm.Blues, \
                    scatter_alpha = 1.0, surface_alpha = 0.6, \
                    gif=False, gif_frames=360, gif_integrals=60):
    ''' v1.1.2  created:2024-03-19  modified:2024-03-24
    (Modified from ---)
    Handle 3D visualization
    '''
    # For full 3D parts
    try:
        #Get max bounds to scale plot effectively; need each axis because of (-50, 50) type bounds
        x_max = xyzarray[:,0].max()
        x_min = xyzarray[:,0].min()
        y_max = xyzarray[:,1].max()
        y_min = xyzarray[:,1].min()
        z_max = xyzarray[:,2].max()
        z_min = xyzarray[:,2].min()
        # Find max total size for each dimension and get the global span
        global_max = max(x_max-x_min, y_max-y_min, z_max-z_min)
        global_pad = global_max * 0.01   # pad dimensions to make plot work better
        if z_max < global_max:
            # Check if z-axis is actually flat-ish; adding a fudge factor in case of future offsetting of flat plates
            #    offsets should still hover around 1 for the average below, so a tight band 0.98-1.02 is fine
            if sum(xyzarray[:,2])/xyzarray.shape[0]>0.98 and sum(xyzarray[:,2])/xyzarray.shape[0]<1.02:
                z_min= z_max  
            # Re-scale z to global
            else:
                z_max = global_max + z_min
                
        # Draw the plots
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(projection='3d')
        ax.scatter3D(xyzarray[:,0], xyzarray[:,1], xyzarray[:,2], s=10, c=point_color, alpha=scatter_alpha)
        ax.plot_trisurf(xyzarray[sphIndices,0], xyzarray[sphIndices,1], xyzarray[sphIndices,2], \
                        triangles=triangIndices, cmap= map_color, alpha=surface_alpha)
        ax.set_xlim3d(x_min - global_pad, x_max + global_pad)
        ax.set_ylim3d(y_min - global_pad, y_max + global_pad)
        ax.set_zlim3d(z_min - global_pad, z_max + global_pad)   #assume Y-Z has already been flipped in xyzarray
        print(f"X dimensions: {x_min},{x_max},{x_max-x_min}")
        print(f"Y dimensions: {y_min},{y_max},{y_max-y_min}")
        print(f"Z dimensions: {z_min},{z_max},{z_max-z_min}")
        
    # for flat plates
    except IndexError:
        #Get max bounds to scale plot effectively; need each axis because of (-50, 50) type bounds
        x_max = xyzarray[:,0].max()
        x_min = xyzarray[:,0].min()
        y_max = xyzarray[:,1].max()
        y_min = xyzarray[:,1].min()
        # Find max total size for each dimension and get the global span
        global_max = max(x_max-x_min, y_max-y_min)
        # if one dimension is bigger, scale the other axis by that amount
        if x_max < global_max:
            # Re-scale z to global
            x_max = global_max + x_min
        if y_max < global_max:
            # Re-scale z to global
            y_max = global_max + y_min
        # pad dimensions to make plot work more good
        global_pad = global_max * 0.01 
            
        # Draw the plots
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(projection='3d')
        ax.scatter3D(xyzarray[:,0], xyzarray[:,1], 0, s=10, c=point_color, alpha=scatter_alpha)
        ax.plot_trisurf(xyzarray[sphIndices,0], xyzarray[sphIndices,1], xyzarray[sphIndices,2], \
                        triangles=triangIndices, cmap= map_color, alpha=surface_alpha)
        ax.set_xlim3d(x_min - global_pad, x_max + global_pad)
        ax.set_ylim3d(y_min - global_pad, y_max + global_pad)
        ax.set_zlim3d(-1, 1)   #assume Y-Z has already been flipped in xyzarray
        print(f"X dimensions: {x_min},{x_max},{x_max-x_min}")
        print(f"Y dimensions: {y_min},{y_max},{y_max-y_min}")
        print(f"Z dimensions: {z_min},{z_max},{z_max-z_min}")
    
        
    # From https://stackoverflow.com/questions/18344934/animate-a-rotating-3d-graph-in-matplotlib:
    def init():
        ax.view_init(elev=10., azim=0)
        return [fig]
        
    def animate(i):
        ax.view_init(elev=10., azim=i)
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
        plt.show()
    

def point_line_distance(point, line_segment_start, line_segment_end):
    '''
    Description: calculate the min distance between a 3D line segment and a point in 3D space; works in 2D as well, but coordinate systems can't be mixed.
    taken almost verbatim from https://stackoverflow.com/questions/56463412/distance-from-a-point-to-a-line-segment-in-3d-python
    '''
    
    #Make sure point and arrays are numpy arrays to allow for operator casting (i.e. '-' operator)
    line_segment_start = np.array(line_segment_start)
    line_segment_end = np.array(line_segment_end)
    point = np.array(point)
    
    # normalized tangent vector
    d = np.divide(line_segment_end - line_segment_start, np.linalg.norm(line_segment_end - line_segment_start))

    # signed parallel distance components
    s = np.dot(line_segment_start - point, d)
    t = np.dot(point - line_segment_end, d)

    # clamped parallel distance
    h = np.maximum.reduce([s, t, 0])

    # perpendicular distance component
    c = np.cross(point - line_segment_start, d)

    return np.hypot(h, np.linalg.norm(c))


def line_line_distance(a0,a1,b0,b1,\
                        clampAll=True, \
                        clampA0=False,clampA1=False,\
                        clampB0=False,clampB1=False):
    ''' v0.2.0  created:2024-05-20  modified:2026-02-11 
    taken verbatim from https://stackoverflow.com/questions/2824478/shortest-distance-between-two-line-segments and modified as needed
        
    INPUT:   Two 'a' and 'b' points; point a/b arbitrary, as is 0/1
                'Clamp' options constrain distances within those segment bound(s); otherwise, the shortest line distance is calculated.
    ACTION:  lorem
    OUTPUT:  
        'pA'        Closest point on segment A
        'pB'        Closest point on segment B
        'distance'  Distance between the points
    '''
    
    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0=True
        clampA1=True
        clampB0=True
        clampB1=True

    #Make sure points are np.arrays
    if type(a0) == list:
        a0 = np.array(a0)
    if type(a1) == list:
        a1 = np.array(a1)
    if type(b0) == list:
        b0 = np.array(b0)
    if type(b1) == list:
        b1 = np.array(b1)
    
    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0
    magA = np.linalg.norm(A)
    magB = np.linalg.norm(B)
        
    _A = A / magA
    _B = B / magB
        
    cross = np.cross(_A, _B)
    denom = np.linalg.norm(cross)**2
        
        
    # If lines are parallel (denom=0) test if lines overlap.
    # If they don't overlap then there is a closest point solution.
    # If they do overlap, there are infinite closest positions, but there is a closest distance
    if not denom:
        d0 = np.dot(_A,(b0-a0))
            
        # Overlap only possible with clamping
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A,(b1-a0))
                
            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 and clampB1:
                    if np.absolute(d0) < np.absolute(d1):
                        return a0,b0,np.linalg.norm(a0-b0)
                    return a0,b1,np.linalg.norm(a0-b1)
                    
                    
            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return a1,b0,np.linalg.norm(a1-b0)
                    return a1,b1,np.linalg.norm(a1-b1)
                    
                    
        # Segments overlap, return distance between parallel segments
        return None,None,np.linalg.norm(((d0*_A)+a0)-b0)
            
        
        
    # Lines criss-cross: Calculate the projected closest points
    t = (b0 - a0);
    detA = np.linalg.det([t, _B, cross])
    detB = np.linalg.det([t, _A, cross])
    
    t0 = detA/denom
    t1 = detB/denom
    
    pA = a0 + (_A * t0) # Projected closest point on segment A
    pB = b0 + (_B * t1) # Projected closest point on segment B
    
    
    # Clamp projections
    if clampA0 or clampA1 or clampB0 or clampB1:
        if clampA0 and t0 < 0:
            pA = a0
        elif clampA1 and t0 > magA:
            pA = a1
            
        if clampB0 and t1 < 0:
            pB = b0
        elif clampB1 and t1 > magB:
            pB = b1
                
        # Clamp projection A
        if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
            dot = np.dot(_B,(pA-b0))
            if clampB0 and dot < 0:
                dot = 0
            elif clampB1 and dot > magB:
                dot = magB
            pB = b0 + (_B * dot)
        
        # Clamp projection B
        if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
            dot = np.dot(_A,(pB-a0))
            if clampA0 and dot < 0:
                dot = 0
            elif clampA1 and dot > magA:
                dot = magA
            pA = a0 + (_A * dot)

    distance = np.linalg.norm(pA-pB)
    
        
    return pA,pB,distance

    
###############################################################################
####   Utility functions   ####################################################
###############################################################################


# Helper function to interpolate between curve points to get a reasonable result
def curve_interpolate(xs, zs, desired_number_of_points,
                        interpol_type = 'quadratic'):
    ''' v1.0   created:2024-03-20   modified:2024-03-20
    Description: Use 2D curve definition to interpolate to arbitrary number of points.
        Obviously breaks down if there are too few points in 'xs' and 'ys', but usually not a problem.
        NOTE: in general for DIW and other additive applications, Z-axis is accurate, but X/Y is arbitrary
            i.e. a 2D slice defining the substrate of a print has a 'real' up/down axis, but many X/Y combos for sideways

    INPUTS:
        'xs'-   in-plane-of-floor coordinates; 'lats', or 'sideways-ness'
        'zs'-   out-of-plane coordinates; 'verts', or 'up/down-ness'
        'desired_number_of_points'- number of interpolated points to return
        (optional)
        'interpol_type'- scipy.interp1d() takes:
                'slinear'- 1st order spline
                'cubic'- 2nd order spline
                'quadratic'- 3rd order spline
    OUTPUT:
        numpy arrays of interpolated xs and zs

    '''

    #TODO: make sure this is a numpy array or list
    x_min = min(xs)
    x_max = max(xs)

    interpol_function = interp1d(xs, zs, kind = interpol_type)

    xsnew = np.linspace(x_min, x_max, desired_number_of_points)
    zsnew = interpol_function(xsnew)

    return xsnew, zsnew


def draw_2D_strand_line_bythickness(interior_points,
                        angle,
                        strand_diam,
                        array_dim,
                        pix_to_um_conv = 1,
                        length = None,
                        line_type = 'simple',
                        thickness_fcn = 'cylinder',
                        show_points = False,
                        show_array_iterations = False,
                        show_final_array = True):
    """
    Description:
        Interpret line points as an ideal strand of radius 'strand_radius' and draw thickness on array with values of microns. 
    INPUTS:
        'interior_points'   array, tuple, or list of iterables; interpret as tuple of array indices (Y,X)
        'strand_diam'       float or int; stand radius in number of pixels (i.e. array indices distance)
                            NOTE: thickness will be 'strand_diam' at center of line; 
        'angle'             float or int; angle between 'interior_point' and array horizontal axis
                            if 'length' is None or 0, assume line is drawn across entire array with center at 'interior_point'
                            otherwise, assume line is drawn in only one direction from 'interior_point' at 'angle'
        'array_dim'         int, tuple, or list; dimensions of 2D array to draw line on
      (optional)
        'pix_to_um_conv'    float or int; side-length of a pixel in microns; used to convert radius      
        'length'            float or int; if NONE, assume line fills the screen edge-to-edge; othewise, assume 'interior_point' is the starting point
                            NOTE: only option for drawing a line in both directions is if 'length'= None. Otherwise 'angle' determines a 2-D line direction.
        'line_type'         str; 'simple'- 2D line from point-to-point
                            'bezier'- bezier curve; control points will be interpolated
        'drop_off_function' str; function defining how 'cylinder' is only function implemented currently
    ACTIONS:
        -lorem
    OUTPUTS:
        'control_point_dict'    dict; specifics of line end points and the 
    """

    #Initialize variables
    control_point_dict = {
        'line_dicts': {},
        'interior_only': True,
        'drawn_array': None
        }
    line_dict = {
        'start_coordinates': [],
        'end_coordinates': [],
        'radius_value': [],
        }
    strand_radius = strand_diam/2

    #TODO: add options in the future
      #if 'strand_thickness_func' is already a callable object, assume it's a function and just use it
    if callable(thickness_fcn):
        pass
    #otherwise, start parsing what the input wants the functional form to look like
      #standard circular strand assumptions
    elif thickness_fcn == 'cylinder':
        strand_thickness_func = lambda distance: 2* math.sqrt(strand_radius**2 - distance**2)
    
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

    #Get array coordinates for drawing the lines and set conditions
      # if no length, go to edges of array from 
    if (not length) and (line_type == 'simple'):
        #In this case, assume a line should be drawn across the entire screen
        control_point_dict['interior_only'] = False
        
        #Run through each passed point
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
            # For right side of point
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

            # For the left side of the point
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

            #Initialize dict for these values and save to output dict
            this_line_dict = line_dict.copy()
            this_line_dict['start_coordinates'] = left_edge_point
            this_line_dict['end_coordinates'] = right_edge_point
            this_line_dict['radius_value'] = strand_radius
            control_point_dict['line_dicts'][f"line_{idx+1}_0"] = this_line_dict

    elif (not length) and (line_type == 'bezier'):
        #TODO: add special bezier flags; not currently implemented
        control_point_dict['interior_only'] = False
    
    elif length and (line_type == 'simple') and (len(interior_points)>=1):

        for idx, start_point in enumerate(interior_points):
            #Find endpoint and coerce to array dimensions
            y_change = start_point[0] * math.sin(math.radians(angle)) *-1  #Flip 'y_change' to match traditional system ((Y,X) with origin at top-left of image)
            x_change = start_point[1] * math.cos(math.radians(angle))
            end_point = [[round(start_point[0] + y_change), 
                            round(start_point[1] + x_change)]]
            #Initialize dict for these values and save to output dict 
            this_line_dict = line_dict.copy()
            this_line_dict['start_coordinates'] = start_point
            this_line_dict['end_coordinates'] = end_point
            this_line_dict['radius_value'] = strand_radius
            control_point_dict['line_dicts'][f"line_{idx+1}_0"] = this_line_dict

      # line starts and ends are now defined
    #Actually draw the thickness onto the return array
    for line_key in list(control_point_dict['line_dicts'].keys()):
        #Pull values to draw
        this_line_dict = control_point_dict['line_dicts'][line_key]
        this_start_point = this_line_dict['start_coordinates']
        this_end_point = this_line_dict['end_coordinates']
        this_radius = this_line_dict['radius_value']

        #Initialize center of line
        if (not length) and (line_type == 'simple'):
            #Draw initial line thickness
            rr, cc = line(int(this_start_point[0]), int(this_start_point[1]), int(this_end_point[0]), int(this_end_point[1]))
            return_array[rr,cc] = this_radius*2
        
            #Find next-index over on each side of centerline, calculate thickness, and add to 

        #TODO: implement bezier for smoother corners
        elif (not length) and (line_type == 'bezier'):
            pass

    return control_point_dict

    


def get_radial_neighbors_array_data(initial_start_coord,
                                   initial_end_coord,
                                   angle,
                                   strand_diameter,
                                   array_dim,
                                   pix_to_um_conv = 1,
                                   thickness_fcn = 'cylinder'):
    """
    Description:
        Handle stepping 1 pixel from centerline to produce a dict of line coordinates and associated thicknesses for radial neighboring points 
        on a strand. I.e. the strand centerline is full thickness, walk 1 pixel to each side and draw a line of that thickness until you reach the
        edge of the strand.

        NOTE: There's probably a better way to do this (at least more elegant), but handling index edge cases is a pain.
    
    INPUTS:
        ''
        ''
        ''
    ACTIONS:
        - 
    OUTPUTS:
        '' 
    """

    #Initialize variables
    starting_points_list = []
    ending_points_list = []
    thickness_list = []
    #Check that start-end points haven't hit an edge-limit
    #  i.e. make sure they aren't continuing on the same edge
    point_left_edge_hit = False
    point_right_edge_hit = False
    #Pull array dimensions
    array_y_dim = array_dim[0]
    array_x_dim = array_dim[1]
    strand_radius = strand_diameter/2
    #Not sure whether 'end' or 'start' is to the right, so find that out and call the 'right' the end
    x_diff = initial_end_coord[1]-initial_start_coord[1]
    y_diff = initial_end_coord[0]-initial_start_coord[0]
    if (x_diff > 0) and (y_diff>=1):
        slope = (initial_end_coord[0]-initial_start_coord[0])/(initial_end_coord[1]-initial_start_coord[1])
        right_coord = initial_end_coord
        left_coord = initial_start_coord
        if slope>0:
            line_orient = 'up'
        elif slope<0:
            line_orient = 'down'
    elif (x_diff < 0) and (y_diff>=1):
        slope = (initial_start_coord[0]-initial_end_coord[0])/(initial_start_coord[1]-initial_end_coord[1])
        right_coord = initial_start_coord
        left_coord = initial_end_coord
        if slope>0:
            line_orient = 'up'
        elif slope<0:
            line_orient = 'down'
      # purely vertical
    elif x_diff == 0:
        slope = 1e8
        line_orient = 'vertical'
        #Convention is that top coord is 'RIGHT'
        if initial_start_coord[0] > initial_end_coord[0]:
            right_coord = initial_start_coord
            left_coord = initial_end_coord
        else:
            right_coord = initial_end_coord
            left_coord = initial_start_coord
      # purely horizontal; already handled the case if 'x_diff==0' above
    elif y_diff <1:
        slope = 0
        line_orient = 'horizontal'
        #Flip coordinates if slope is negative
        if x_diff>0:
            right_coord = initial_end_coord
            left_coord = initial_start_coord
        elif x_diff<0:
            right_coord = initial_start_coord
            left_coord = initial_end_coord

    #Get original line equation


    #TODO: add options in the future for other thickness drop-off shapes 
      #if 'strand_thickness_func' is already a callable object, assume it's a function and just use it
    if callable(thickness_fcn):
        pass
    #otherwise, start parsing what the input wants the functional form to look like
      #standard circular strand assumptions
    elif thickness_fcn == 'cylinder':
        strand_thickness_func = lambda distance: 2* math.sqrt(strand_radius**2 - distance**2)

    #####################################################
    #Check for edge coordinates
      # LEFT point/start coordinates
    left_y = left_coord[0]
    left_x = left_coord[1]
    if (left_y == 0) and ((left_x>0) and (left_x < array_x_dim)):
        left_edge_pos = 'north'
    elif (left_y == 0) and (left_x == 0):
        left_edge_pos = 'northwest'
    elif (left_y == 0) and (left_x== (array_x_dim-1)):
        left_edge_pos = 'northeast'

    elif (left_y == (array_y_dim-1)) and ((left_x>0) and (left_x < array_x_dim)):
        left_edge_pos = 'south'
    elif (left_y == (array_y_dim-1)) and (left_x == 0):
        left_edge_pos = 'southwest'
    elif (left_y == (array_y_dim-1)) and (left_x== (array_x_dim-1)):
        left_edge_pos = 'southeast'

    elif (left_x == 0) and ((left_y>0) and (left_y < array_y_dim)):
        left_edge_pos = 'west'
    elif (left_x == 0) and (left_y == 0):
        left_edge_pos = 'northwest'
    elif (left_x == 0) and (left_y== (array_y_dim-1)):
        left_edge_pos = 'southwest'

    elif (left_x == (array_x_dim-1)) and ((left_y>0) and (left_y < array_y_dim)):
        left_edge_pos = 'east'
    elif (left_x == (array_x_dim-1)) and (left_y == 0):
        left_edge_pos = 'northeast'
    elif (left_x == (array_x_dim-1)) and (left_y== (array_y_dim-1)):
        left_edge_pos = 'southeast'
    else:
        left_edge_pos = 'interior'

      # RIGHT point/end coordinates
    right_y = right_coord[0]
    right_x = right_coord[1]
    if (right_y == 0) and ((right_x>0) and (right_x < array_x_dim)):
        right_edge_pos = 'north'
    elif (right_y == 0) and (right_x == 0):
        right_edge_pos = 'northwest'
    elif (right_y == 0) and (right_x== (array_x_dim-1)):
        right_edge_pos = 'northeast'

    elif (right_y == (array_y_dim-1)) and ((right_x>0) and (right_x < array_x_dim)):
        right_edge_pos = 'south'
    elif (right_y == (array_y_dim-1)) and (right_x == 0):
        right_edge_pos = 'southwest'
    elif (right_y == (array_y_dim-1)) and (right_x== (array_x_dim-1)):
        right_edge_pos = 'southeast'

    elif (left_x == 0) and ((right_y>0) and (right_y < array_y_dim)):
        right_edge_pos = 'west'
    elif (left_x == 0) and (right_y == 0):
        right_edge_pos = 'northwest'
    elif (left_x == 0) and (right_y== (array_y_dim-1)):
        right_edge_pos = 'southwest'

    elif (left_x == (array_x_dim-1)) and ((right_y>0) and (right_y < array_y_dim)):
        right_edge_pos = 'east'
    elif (left_x == (array_x_dim-1)) and (right_y == 0):
        right_edge_pos = 'northeast'
    elif (left_x == (array_x_dim-1)) and (right_y== (array_y_dim-1)):
        right_edge_pos = 'southeast'
    else:
        right_edge_pos = 'interior'

    #####################################################

    #If points are on edges, proceeed with edge calcs
    if (right_edge_pos != 'interior') and (left_edge_pos != 'interior'):
        radial_distance = 0
        last_right_right = right_coord
        last_right_left = right_coord
        last_left_right = left_coord
        last_left_left = left_coord

        #At each step: 
        #   1) move 1 pixel in each direction from strand center
        #   2) adjust for any crossing of corners or (for interior points) running into edges,
        #   3) calculate an 'average' radial distance for both points
        #   4) add points and distance to lists for return
        #   5) increment distance until a line is drawn that is beyond the strand radius
        iteration_cntr = 0
        max_iterations = int(max(array_x_dim, array_y_dim))
        while (radial_distance < strand_radius) and (iteration_cntr<max_iterations):
            #NOTE: 'Left' incremented down in X and down in Y; 'Right' increment
            #
            #      Confusing names here are because there are two points (left and right) defining the initial line.
            #      In order to find next-line-over, we need to increment both points with new points on either side of the initial line.
            #      This leads 'right' and 'left' center points to have two new points that are each propogated left and right.
            #      First split is above.
            
            #If top or bottom, shift 1 pixel in each direction and calculate effect
            if (right_edge_pos == 'north') or (right_edge_pos == 'south'):
                
                #If slope is positive, things are easy
                if slope > 0:
                    #RIGHT-right
                      # try basic increment right
                    if (last_right_right[1]+1 >= 0) and (last_right_right[1]+1 < array_x_dim):
                        new_right_right_coord = [last_right_right[0], last_right_right[1]+1]
                      # running off the edge north/south-east; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_right_right[1]+1 >= array_x_dim) and (slope > 1e7):
                        new_right_right_coord = last_right_right
                      # running off the edge north/south-east; move to 'east' edge
                    elif (last_right_right[1]+1 >= array_x_dim):
                        new_right_right_coord = [last_right_right[0]+1, last_right_right[1]]
                        right_edge_pos = 'east'
                    
                    #RIGHT-left
                      # try basic increment left
                    if (last_right_left[1]-1 >= 0) and (last_right_left[1]-1 < array_x_dim):
                        new_right_left_coord = [last_right_left[0], last_right_left[1]-1]
                      # running off the edge north/south-west; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_right_left[1]-1 <0) and (slope > 1e7):
                        new_right_left_coord = last_right_left
                      # running off the edge north/south-west; (+) slope, so assume there's no other place to go and keep the last position
                    elif (last_right_left[1]-1 <0):
                        new_right_left_coord = last_right_left
                
                #If slope is negative, 'X' and 'Y' increments move in different directions
                #  i.e. 'left' shift is +X and 'right' shift is -X  
                if slope < 0:
                    #RIGHT-right
                      # try basic increment right
                    if (last_right_right[1]-1 >= 0) and (last_right_right[1]-1 < array_x_dim):
                        new_right_right_coord = [last_right_right[0], last_right_right[1]+1]
                      # running off the edge north/south-west; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_right_right[1]-1 < 0) and (slope > 1e7):
                        new_right_right_coord = last_right_right
                      # running off the edge north/south-west; (-) slope, so assume there's no other place to go and keep the last position
                    elif (last_right_right[1]-1 < 0):
                        new_right_right_coord = last_right_right
                    
                    #RIGHT-left
                      # try basic increment left
                    if (last_right_left[1]+1 >= 0) and (last_right_left[1]+1 < array_x_dim):
                        new_right_left_coord = [last_right_left[0], last_right_left[1]-1]
                      # running off the edge north/south-east; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_right_left[1]+1 >= array_x_dim) and (slope > 1e7):
                        new_right_left_coord = last_right_left
                      # running off the edge north/south-east; shift to 'east' position and increment -Y instead of +X
                    elif (last_right_left[1]+1 >= array_x_dim):
                        new_right_left_coord = [last_right_left[0]+1, last_right_left[1]]
                        right_edge_pos = 'east'
                
                #If horizontal, easy case because point shouldn't be ; convention is 'left' shift is UP (-Y)
                if slope == 0:
                    pass

            if (left_edge_pos == 'north') or (left_edge_pos == 'south'):
                #If slope is positive, things are easy
                if slope > 0:
                    #LEFT-right
                      # try basic increment right
                    if (last_left_right[1]+1 >= 0) and (last_left_right[1]+1 < array_x_dim):
                        new_left_right_coord = [last_left_right[0], last_left_right[1]+1]
                      # running off the edge north/south-east; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_left_right[1]+1 >= array_x_dim) and (slope > 1e7):
                        new_left_right_coord = last_left_right
                      # running off the edge north/south-east; (+) slope, so assume there's no other place to go and keep the last position
                    elif (last_left_right[1]+1 >= array_x_dim):
                        new_left_right_coord = last_left_right
                    
                    #LEFT-left
                      # try basic increment left
                    if (last_left_left[1]-1 >= 0) and (last_left_left[1]-1 < array_x_dim):
                        new_left_left_coord = [last_left_left[0], last_left_left[1]-1]
                      # running off the edge north/south-west; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_left_left[1]-1 <0) and (slope > 1e7):
                        new_left_left_coord = last_left_left
                      # running off the edge north/south-west; (+) slope, so assume there's no other place to go and keep the last position
                    elif (last_left_left[1]-1 <0):
                        new_left_left_coord = last_left_left
                
                #If slope is negative, 'X' and 'Y' increments move in different directions
                #  i.e. 'left' shift is +X and 'right' shift is -X  
                if slope < 0:
                    #LEFT-right
                      # try basic increment right
                    if (last_left_right[1]-1 >= 0) and (last_left_right[1]-1 < array_x_dim):
                        new_left_right_coord = [last_left_right[0], last_left_right[1]-1]
                      # running off the edge north/south-west; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_left_right[1]-1 < 0) and (slope > 1e7):
                        new_left_right_coord = last_left_right
                      # running off the edge north/south-west; shift to 'west' position and start walking down the 'west' edge
                    elif (last_left_right[1]-1 < 0):
                        new_left_right_coord = [last_left_right[0]-1, last_left_right[1]]
                        right_edge_pos = 'west'
                    
                    #LEFT-left
                      # try basic increment left
                    if (last_left_left[1]+1 >= 0) and (last_left_left[1]+1 < array_x_dim):
                        new_left_left_coord = [last_left_left[0], last_left_left[1]-1]
                      # running off the edge north/south-east; vertical line, so just pass the last line again (nowhere to increment to)
                    elif (last_left_left[1]+1 >= array_x_dim) and (slope > 1e7):
                        new_left_left_coord = last_left_left
                      # running off the edge north/south-east; (-) slope, so assume there's nowhere else to go
                    elif (last_left_left[1]+1 >= array_x_dim):
                        new_left_left_coord = last_left_left
                
                #If horizontal, easy case because case shouldn't exist ; convention is "'left' shift is UP (-Y)"
                if slope == 0:
                    pass

            #If left or right, shift 1 pixel in each direction and calculate effect
            if (right_edge_pos == 'east') or (right_edge_pos == 'west'):
                #If slope is positive, things are easy
                if slope > 0:
                    #RIGHT-right
                      # try basic increment right
                    if (last_right_right[0]+1 >= 0) and (last_right_right[0]+1 < array_y_dim):
                        new_right_right_coord = [last_right_right[0]+1, last_right_right[1]]
                      # running off the edge south-east; assume there's no where else to go
                    elif (last_right_right[0]+1 >= array_y_dim):
                        new_right_right_coord = last_right_right
                    
                    #RIGHT-left
                      # try basic increment left
                    if (last_right_left[0]-1 >= 0) and (last_right_left[0]-1 < array_y_dim):
                        new_right_left_coord = [last_right_left[0]-1, last_right_left[1]]
                      # running off the north-east edge; transition to north
                    elif (last_right_left[0]-1 <0):
                        new_right_left_coord = [last_right_left[0], last_right_left[1]-1]
                        right_edge_pos = 'north'

                #If slope is negative, 'X' and 'Y' increments move in different directions
                #  i.e. 'left' shift is +X and 'right' shift is -X  
                if slope < 0:
                    #RIGHT-right
                      # try basic increment right
                    if (last_right_right[0]+1 >= 0) and (last_right_right[0]+1 < array_y_dim):
                        new_right_right_coord = [last_right_right[0], last_right_right[1]+1]
                      # running off the edge south-west; no where else to go, just pass last value
                    elif (last_right_right[0]+1 >= array_y_dim):
                        new_right_right_coord = last_right_right
                    
                    #RIGHT-left
                      # try basic increment left
                    if (last_right_left[0]-1 >= 0) and (last_right_left[0]-1 < array_y_dim):
                        new_right_left_coord = [last_right_left[0], last_right_left[1]-1]
                      # running off the edge north/south-east; shift to 'east' position and increment -Y instead of +X
                    elif (last_right_left[0]-1 <0):
                        new_right_left_coord = [last_right_left[0], last_right_left[1]-1]
                        right_edge_pos = 'north'
                
                #If horizontal, easy case; convention is 'left' shift is UP (-Y)
                if slope == 0:
                    #RIGHT-right
                    if (last_right_right[0]+1 >= 0) and (last_right_right[0]+1 < array_y_dim):
                        new_right_right_coord = [last_right_right[0], last_right_right[1]+1]
                    elif (last_right_right[0]+1 >= array_y_dim):
                        new_right_right_coord = last_right_right
                    #RIGHT-left
                    if (last_right_left[0]-1 >= 0) and (last_right_left[0]-1 < array_y_dim):
                        new_right_left_coord = [last_right_left[0], last_right_left[1]-1]
                    elif (last_right_left[0]-1 < 0):
                        new_right_left_coord = last_right_left

            #If left or right, shift 1 pixel in each direction and calculate effect
            if (left_edge_pos == 'east') or (left_edge_pos == 'west'):
                #If slope is positive, things are easy
                if slope > 0:
                    #LEFT-right
                      # try basic increment right
                    if (last_left_right[0]+1 >= 0) and (last_left_right[0]+1 < array_y_dim):
                        new_left_right_coord = [last_left_right[0]+1, last_left_right[1]]
                      # running off the edge south-west; transition to 'south' 
                    elif (last_left_right[0]+1 >= array_y_dim):
                        new_left_right_coord = [last_left_right[0], last_left_right[1]+1]
                        left_edge_pos = 'south'
                    
                    #LEFT-left
                      # try basic increment left
                    if (last_left_left[0]-1 >= 0) and (last_left_left[0]-1 < array_y_dim):
                        new_left_left_coord = [last_left_left[0]-1, last_left_left[1]]
                      # running off the north-east edge; transition to north
                    elif (last_left_left[0]-1 <0):
                        new_left_left_coord = [last_left_left[0], last_left_left[1]+1]
                        right_edge_pos = 'north'

                #If slope is negative, 'X' and 'Y' increments move in different directions
                #  i.e. 'left' shift is +X and 'right' shift is -X  
                if slope < 0:
                    #LEFT-right
                      # try basic increment right
                    if (last_left_right[0]+1 >= 0) and (last_left_right[0]+1 < array_y_dim):
                        new_left_right_coord = [last_left_right[0]+1, last_left_right[1]]
                      # running off the edge south-west; no where else to go, just pass last value
                    elif (last_left_right[0]+1 <0):
                        new_left_right_coord = last_left_right
                    
                    #LEFT-left
                      # try basic increment left
                    if (last_left_left[0]-1 >= 0) and (last_left_left[0]-1 < array_y_dim):
                        new_left_left_coord = [last_left_left[0]-1, last_left_left[1]]
                      # running off the edge north/south-east; shift to 'east' position and increment -Y instead of +X
                    elif (last_left_left[0]-1 <0):
                        new_left_left_coord = [last_left_left[0]-1, last_left_left[1]]
                        right_edge_pos = 'north'
                
                #If horizontal, easy case; convention is 'left' shift is UP (-Y)
                if slope == 0:
                    #LEFT-right
                    if (last_left_right[0]+1 >= 0) and (last_left_right[0]+1 < array_y_dim):
                        new_left_right_coord = [last_left_right[0], last_left_right[1]+1]
                    elif (last_left_right[0]+1 >= array_y_dim):
                        new_left_right_coord = last_left_right
                    #LEFT-left
                    if (last_left_left[0]-1 >= 0) and (last_left_left[0]-1 < array_y_dim):
                        new_left_left_coord = [last_left_left[0], last_left_left[1]-1]
                    elif (last_left_left[0]-1 < 0):
                        new_left_left_coord = last_left_left

            #Calculate distance between initial centerline and this new line
            left_distances = [point_line_distance(new_left_left_coord, initial_start_coord, initial_end_coord),
                              point_line_distance(new_right_left_coord, initial_start_coord, initial_end_coord)]
            average_left_distance = sum(left_distances)/2
            right_distances = [point_line_distance(new_left_right_coord, initial_start_coord, initial_end_coord),
                              point_line_distance(new_right_right_coord, initial_start_coord, initial_end_coord)]
            average_left_distance = sum(left_distances)/2
            average_right_distance = sum(right_distances)/2

            radial_distance = max(average_left_distance,average_right_distance) * pix_to_um_conv            
                        
            #Store values if progression isn't continuing down the same array edge (store first instance and no others)
            #NOTE: order doesn't matter here, but trying to maintain 'left = start', 'right = end' convention for (some) clarity
              #new ->right line
            if not point_right_edge_hit:
                starting_points_list.append(new_left_right_coord)
                ending_points_list.append(new_right_right_coord)
                thickness_list.append(average_right_distance)
            if not point_left_edge_hit:
                starting_points_list.append(new_left_left_coord)
                ending_points_list.append(new_right_left_coord)
                thickness_list.append(average_left_distance)
            
            #Check to make sure a line isn't being drawn along the edge again (allowed first time)
            left_edge_bool = bool(
                ((new_left_left_coord[0] == 0) and (new_right_left_coord[0] ==0)) or
                ((new_left_left_coord[0] == array_y_dim) and (new_right_left_coord[0] == array_y_dim)) or
                ((new_left_left_coord[1] == 0) and (new_right_left_coord[1] == 0)) or
                ((new_left_left_coord[1] == array_x_dim) and (new_right_left_coord[1] == array_x_dim))
                )
            right_edge_bool = bool(
                ((new_left_right_coord[0] == 0) and (new_right_right_coord[0] ==0)) or
                ((new_left_right_coord[0] == array_y_dim) and (new_right_right_coord[0] == array_y_dim)) or
                ((new_left_right_coord[1] == 0) and (new_right_right_coord[1] == 0)) or
                ((new_left_right_coord[1] == array_x_dim) and (new_right_right_coord[1] == array_x_dim))
                )
              # change edge flag if just hit
            if (left_edge_bool) and (not point_left_edge_hit):
                point_left_edge_hit = True
            if (right_edge_bool) and (not point_right_edge_hit):
                point_right_edge_hit = True

            #Store current position for next iteration
            last_right_right = new_right_right_coord
            last_right_left = new_right_left_coord
            last_left_right = new_left_right_coord
            last_left_left = new_left_left_coord

            #Increment cntr in case things go poorly above
            iteration_cntr += 1


    #If points are interior, proceed with interior calcs
    elif (right_edge_pos == 'interior') and (left_edge_pos == 'interior'):
        pass
    
    #If points are mixed, do both interior and edge calcs
    #TODO: fill in these cases
      # end pos only edge
    elif (right_edge_pos != 'interior'):
        pass
      # start pos only edge
    elif (left_edge_pos != 'interior'):
        pass


    neighbor_dict = {
        'starting_points': starting_points_list,
        'ending_points': ending_points_list,
        'thicnkess': thickness_list
        }

    return neighbor_dict