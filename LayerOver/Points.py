# -*- coding: utf-8 -*-
"""
© 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security 
Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of 
Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf 
a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

Created:   2024-03-18
Modified:  2025-07-15
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

    

