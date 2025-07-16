# -*- coding: utf-8 -*-
"""
Created:   2025-07-16
Modified:  2025-07-16
Version:   0.1.0

@author: Aaron Pital (Los Alamos National Lab)

Description: 
    
    
Change log:
    2025-07-16 (initial commit)
        - Split off POV and other voxel rendering to streamline import namespace.

"""

class Camera():

    def __init__(self, gcode_part, lens_dict = None, render_dict= None, view_direction = 'in'):
        ''' v0.2.0   created:2024-04-10   modified:2024-05-08

        INPUTS:  'cart_location_dict' - dictionary of points to view from; need location and normal vector
                 'render_dict' - contains all the geometry information
                 'strand_geometry' - 'by-layer'-- each layer has it's own strand geometry
        '''

        self.gcode_part = gcode_part
        
        # If 'lens_dict' isn't provided, make one
        if lens_dict == None:
            self.lens_dict = {
                'type': 'simple',
                'shape': 'square',
                'profile': 'flat',
                'size': (6,6),
                'working_distance': 6
                }
        else:
            self.lens_dict = lens_dict

          #Modify with gcode_part
        

        # If 'render_dict' isn't specified, make one
        if render_dict == None:
            self.render_dict = {
                'geometry_dict':
                    {
                    'strand_size_type': 'global',
                    'strand_diameter': 150,
                    },
                
                'intensity_cmap': None,
                'perspective_padding': 1
                }

        #Instantiate view dictionaries
        self.povs = {}
        self.voxels = {}

    
    def make_povs(self, number_of_locations, location_list =None):
        ''' v0.1.0   created:2024-05-09   modified:2024-05-20

        INPUTS:  'location_number' - number of points to pull cameras on
                 (optional)
                 'location_list'- supply a list of discrete locations rather than pick at random
        '''

        if location_list == None:
        
            point_dicts = Camera.get_random_point_summaries(self.gcode_part, n_retrieval_points = number_of_locations, \
                                                            excluded_point_indcs = [], show_iterations = True)
            self.pov_point_summaries = point_dicts
        else:
            #TODO: make this part work
            pass

        for idx, key in enumerate(list(point_dicts.keys())):
            center = point_dicts[idx]['approximate_center']
            normal_vector = point_dicts[idx]['unit_normal']
            radius = self.lens_dict['size'][0]
            
            #get voxel bounds and add to dictionary 
            voxel_dict = Camera.get_voxel(center, normal_vector, radius= radius, voxel_depth= 2*radius, n_pixels=1000)
            voxel_dict = Camera.get_layer_points(self, voxel_dict)

            #Updtate 'voxel' dictionary on object
            self.voxels.update( {key:voxel_dict} )
    
    
    def populate_povs(self, make_opacity = True, direction = 'in'):
        
        ''' v0.1.1   created:2024-05-09   modified:2024-05-11

        INPUTS:  'cart_location_dict' - dictionary of points to view from; need location and normal vector
                 'render_dict' - contains all the geometry information
                 'strand_geometry' - 'by-layer'-- each layer has it's own strand geometry
        '''
        
        #Initialize the biz
        point_id_list = list(self.location_dict.keys())
        
        #Run through each point's properties in 'cart_location_dict', get properties and apply to 'self.povs' dictionary
        for point_id in point_id_list:
            this_center = self.location_dict[point_id]['approximate_center']
            this_normal = self.location_dict[point_id]['unit_normal']
            #Derive eqn of plane (cx*nx + cy*ny + cz*nz = d) by plugging in center to normal
            plane_d = this_center[0]*this_normal[0] + this_center[1]*this_normal[1] + this_center[2]*this_normal[2]
              
              #Create
            if direction == 'in':
                pass
            elif direction == 'out':
                pass
            else:
                pass    #default to 
            
            
            if self.render_dict['geometry_dict']['strand_size_type'] == 'global':
                pass
            
            
    @staticmethod
    def make_camera(center_point, unit_normal, 
                    plane_size = (6,6), lens_type = 'square', 
                    focal_pyramid_length = 6, 
                    camera_color= 'k', opacity = 0.7 
                    ):
        ''' v0.1.0  created:2024-05-10  modified:2024-05-10

        Take 
        
            INPUT:   'center'- point at center of camera/lens plane
                     'unit_normal'- unit normal vector for direction of camera
            ACTION:  math
            OUTPUT:  'camera_dict'-
                        'plane_points'- 4 pts (x,y,z) describing render plane
                        'lens_type'- defaults to a square lens 6mm x 6mm
                        'voxel_points'- 8 pts (x,y,z) describing bounding box for camera's observation voxel
                             (includes padding so geometries behave)
                        'focal_pyramid_point'- point 'behind' camera lens that defines focal cone (shown as pyramid)
        '''
        
        #plane
        
        #return camera_dict
        
        pass

