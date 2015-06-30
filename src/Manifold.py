import numpy as np
import pdb
from sklearn.neighbors import NearestNeighbors
#pdb.set_trace()

class Manifold(object):
    """Implments the timestepping
    """
    
    def __init__(self,max_points,number_of_bands):
        self.max_points = max_points
        self.points = np.zeros((max_points,3))
        self.colors = np.zeros((max_points))
        self.triangles = np.zeros((max_points,3),dtype=np.int32)
        self.closest_points = np.zeros((max_points))
        self.number_of_bands = number_of_bands
        self.max_ring = number_of_bands +1
        self.ring_starts = np.zeros((self.max_ring),dtype=np.int32)
        self.current_row = 0
        self.current_ring = 0
        self.ring_starts[0] = 0
        self.current_triangle = 0

    def clear_visual_data(self):
        max_points = self.max_points
        self.colors = np.zeros((max_points))
        self.triangles = np.zeros((max_points,3),dtype=np.int32)
        self.current_triangle = 0
    
    def add_ring(self,band):
        n_rows_in_band = band.shape[0]
        if self.current_ring > 0:
            last_ring_start = self.ring_starts[self.current_ring - 1]
            last_band = self.points[last_ring_start:self.current_row,:]
            nbrs = NearestNeighbors(n_neighbors=1).fit(band)
            distances, indices = nbrs.kneighbors(last_band)
            pdb.set_trace()
            self.closest_points[last_ring_start:self.current_row] = indices[:,0] + self.current_row
        self.points[(self.current_row):(self.current_row + n_rows_in_band),:] = band.copy()
        self.current_row += n_rows_in_band
        self.current_ring += 1
        self.ring_starts[self.current_ring] = self.current_row
    
    def set_colors(self,delete_bands=0):
        ring_num = 0
        while ring_num < self.current_ring - delete_bands:
            begin = self.ring_starts[ring_num]
            end = self.ring_starts[ring_num+1]
            self.colors[begin:end].fill(ring_num)
            ring_num += 1
    
    def add_triangle(self,triangle):
        self.triangles[self.current_triangle,:] = np.array(triangle)
        self.current_triangle += 1
        
    def assign_bottom_triangles(self,delete_bands):
        # Make the bottom triangles in a vectorized manner. These
        # Triangles are of the form
        # [point_i_on_bottom,
        #   point_i+1_on_bottom,
        #   point_closest_to_point_i_thats_on_top_ring]
        first_point_next_ring = self.ring_starts[self.current_ring-1-delete_bands]
        number_of_bottom_triangles = first_point_next_ring
        start_tri = self.current_triangle
        end_tri = self.current_triangle + number_of_bottom_triangles
        self.triangles[start_tri:end_tri,0] = np.arange(number_of_bottom_triangles)
        self.triangles[start_tri:end_tri,1] = np.arange(number_of_bottom_triangles) + 1
        # Change end points to beginning points
        self.triangles[self.ring_starts[1:self.current_ring] - 1,1] = \
            self.ring_starts[0:self.current_ring-1]
        self.triangles[start_tri:end_tri,2] = self.closest_points[0:first_point_next_ring]
        self.current_triangle += number_of_bottom_triangles
    
    def assign_top_triangles(self,delete_bands):
        # Make the bottom triangles in a vectorized manner. These
        # Triangles are of the form
        # [point_i_on_bottom,
        #   point_i+1_on_bottom,
        #   point_closest_to_point_i_thats_on_top_ring]
        first_point_next_ring = self.ring_starts[self.current_ring-1-delete_bands]
        number_of_top_triangles = first_point_next_ring
        start_tri = self.current_triangle
        end_tri = self.current_triangle + number_of_top_triangles
        self.triangles[start_tri:end_tri,1] = np.arange(number_of_top_triangles) + 1
        self.triangles[self.current_triangle + self.ring_starts[1:self.current_ring] - 1,1] = \
            self.ring_starts[0:self.current_ring-1]
        self.triangles[start_tri:end_tri,0] = self.closest_points[0:first_point_next_ring]
        self.triangles[start_tri:end_tri,2] = self.closest_points[self.triangles[start_tri:end_tri,1]]
        self.current_triangle += number_of_top_triangles
        
    #def assign_top_triangles(self):
        # First, we get a list of the points that were involved in triangles
        # that were as the point
        # point_closest_to_point_i_thats_on_top_ring
        # or
        # point_i+1_on_bottom
        # in a triangle of the type
        # [point_i_on_bottom,
        #   point_i+1_on_bottom,
        #   point_closest_to_point_i_thats_on_top_ring]
        #u, indices = np.unique(a, return_index=True)
        # The below is typically what we will want
        #top_points_normal_order, indices_normal_order = \
        #    np.unique(self.triangles[0:self.current_triangle,2], 
        #              return_index=True)
        # Now we get a list of 
        # point_closest_to_point_i_thats_on_top_ring
        # that are associated with 
        # point_i_on_bottom
        # as teh first point in the triangle
        #top_points_from_first_points = self.closest_points[self.ring_starts[0:self.current_ring]]
        # Need to figure out this code
        #top_points_reverse_order, indices_reverse_order = \
        #    np.unique(self.triangles[0:self.current_triangle,2], 
        #              return_index=True)
        # Check for intersection in top_points_from_first_points, top_points_reverse_order
        # Check for difference in first_points, if it's greater than 8, flag it
        
        
        
        
        

    def make_triangles(self,delete_bands=0):
        # Must fix for top ring breaking periodic in the middle!!!!!
        # looks like I should make an iterator for range(closest,next_closest)
        # Should be checking appropraite bounds in the upper ring and working 
        # modulo possibly. At least do some manual trick.
        #
        self.assign_bottom_triangles(delete_bands)
        self.assign_top_triangles(delete_bands)
            
        #point = int(self.ring_starts[ring_num])
        #while point < self.ring_starts[ring_num + 1] - 1:
            #closest = int(self.closest_points[point])
            #self.add_triangle([point,point+1,closest])
            #next_closest = int(self.closest_points[point+1])
            #while (closest == next_closest):
                #point += 1
                #closest = int(self.closest_points[point])
                #next_closest = int(self.closest_points[point+1])
                #self.add_triangle([point,point+1,closest])
            #for this_top in range(closest,next_closest):
                #self.add_triangle([point+1,this_top,this_top+1])
            #point += 1
        #pdb.set_trace()
        ## Now we are on the last point in the ring. The next point on
        ## this ring is also the first point.
        #closest = int(self.closest_points[point])
        #pointplusone = self.ring_starts[ring_num]
        #self.add_triangle([point,pointplusone,closest])
        #next_closest = int(self.closest_points[pointplusone])
        #if closest != next_closest:
            #for this_top in range(closest,self.ring_starts[ring_num + 2] - 2):
                #self.add_triangle([pointplusone, this_top, this_top+1])
            #self.add_triangle([pointplusone, self.ring_starts[ring_num + 2] - 1, self.ring_starts[ring_num + 1]])
            
            
   



#indices                                           
#distances     
    #def __init__(self, bottom_points, top_points):
        #"""Creates a blank ring object using an options object.
        
        #Args:
            #options: An ring options object.
        #"""
        #points_all = np.vstack((bottom_points,top_points))
        #vertex_start = range(bottom_points.shape[0])
        #vertex_end = range(points_all.shape[0]-1,bottom_points.shape[0]-1,-1)
        #first_cell = []
        #first_cell.append(vertex_end[-1])
        #first_cell.extend(vertex_start)
        #first_cell.extend(vertex_end[:-1])
        #tvtk.PolyData.__init__(self,points=points_all,polys=np.array([first_cell]))
        #last_points=top_points
        
        
    #def add_points(self, next_points):
        #"""Sets the first ring as a circle
        
        #"""
        #points_all = np.vstack((self.last_points,next_points))
        #vertex_start = range(self.last_points.shape[0])
        #vertex_end = range(points_all.shape[0]-1,self.last_points.shape[0]-1,-1)
        #this_cell = []
        #this_cell.extend(vertex_start)
        #this_cell.extend(vertex_end)
        #this_cell.append(vertex_start[0])
        #new_obj = tvtk.PolyData(points=points_all,polys=np.array(this_cell))
        #self.copy_cells(new_obj,[0])
        #self.last_points = next_points.copy()
        
