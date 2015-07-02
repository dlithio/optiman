import numpy as np
from sklearn.neighbors import NearestNeighbors
from mayavi import mlab
import inspect,os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
dname += os.sep + 'results'
os.chdir(dname)

class Manifold(object):
    """Implments the timestepping
    """
    
    def __init__(self,folder_name):
        # Load header for loading real file
        header = np.loadtxt(folder_name + '/header')
        ndim = int(header[0])
        # Load an actual rings file
        rings = self.read_array(folder_name +'/output',np.float64)
        numrows = rings.shape[0]/(ndim+1)
        rings.shape = (numrows,ndim+1)
        # Function to get start points of each ring

        start_points = self.get_start_points(rings)
        end_points = []
        for point in start_points[1:]:
            end_points.append(point - 1)
        end_points.append(rings.shape[0] - 1)
        # Make an array that is for the bottom triangles and is (numpoints-points_in_top_ring x 3)
        # Fill in first two columns based on start points of each ring
        bottom_triangles = np.zeros((start_points[-1],3),dtype=np.int)
        bottom_triangles[:,0] = np.arange(0,start_points[-1])
        bottom_triangles[:,1] = np.arange(1,start_points[-1]+1)
        bottom_triangles[end_points[:-1],1] = start_points[:-1]
        # Make an array for the top triangles that is (numpoints-points_in_bottom_ring x 3)
        # Fill in first two columns based on start points of each ring
        top_triangles = np.zeros((rings.shape[0] - start_points[1],3),dtype=np.int)
        top_triangles[:,0] = np.arange(start_points[1],rings.shape[0])
        top_triangles[:,1] = np.arange(start_points[1]+1,rings.shape[0]+1)
        for start_point,end_point in zip(start_points[1:],end_points[1:]):
            top_triangles[end_point-start_points[1],1] = start_point

        # Use the function to build the bottom triangles
        for ring_num in range(len(start_points) - 1):
            bottom_ring_start = start_points[ring_num]
            bottom_ring_end = end_points[ring_num] + 1
            top_ring_start = start_points[ring_num + 1]
            top_ring_end = end_points[ring_num + 1] + 1
            bottom_ring = rings[bottom_ring_start:bottom_ring_end,1:]
            top_ring = rings[top_ring_start:top_ring_end,1:]
            bottom_triangles[bottom_ring_start:bottom_ring_end,2] = \
                         self.get_closest_points(bottom_ring,top_ring,top_ring_start)
                         
        self.get_top_to_bottom(bottom_triangles,top_triangles,start_points,end_points)
        
        all_triangles = np.concatenate((bottom_triangles,top_triangles), axis=0)
        self.points = rings[:,1:]
        self.triangles = all_triangles
        self.time = rings[:,0]
        
        # Now visualize it
        # dims = [1,2,3]
        #mesh = tvtk.PolyData(points=rings[:,dims], polys=all_triangles)
        #self.points = rings[:,dims]
        #mesh.point_data.scalars = rings[:,0]
        #mesh.point_data.scalars.name = 'time'
        #fdot = read_array('fdot',np.float64)
        #fdot.shape = (fdot.shape[0]/2,2)
        #mesh.point_data.scalars = fdot[:,1]
        #mesh.point_data.scalars.name = 'f_dot_t'
        #t_angle = read_array('t_angle',np.float64)
        #t_angle.shape = (t_angle.shape[0]/2,2)
        #mesh.point_data.scalars = t_angle[:,1]
        #mesh.point_data.scalars.name = 't_angle'
        #self.view(mesh)
    
    def draw_manifold(self,dims):
        mlab.triangular_mesh(self.points[:,0], self.points[:,1], self.points[:,2], self.triangles, scalars=self.time)
    
    # Now the test of the new visualization software
    def read_array(self,filename,datatype):
        try:
            result = np.fromfile(filename,dtype=datatype)
        except IOError:
            print(filename + " was not yet available")
        return result
        
    def get_start_points(self,rings):
        start_points = []
        ring_num = 0
        for row_num,row_ring_num in enumerate(rings[:,0]):
            if row_ring_num != ring_num:
                ring_num = row_ring_num
                start_points.append(row_num)
        return start_points
        
    # Next, we want a function that takes (bottom ring, top ring, closest points to bottom on the top)
    #   and returns an array as long as the top ring array that is the point on the bottom ring that we
    #   should connect to when the corresponding point in the top is first in the triangle
    def get_top_to_bottom(self,bottom_triangles,top_triangles,start_points,end_points):
        top_triangles[:,2] = -1
        for triangle in bottom_triangles:
            if triangle[1] > top_triangles[triangle[2] - start_points[1],2]:
                top_triangles[triangle[2] - start_points[1],2] = triangle[1]
        # Make corrections for the first point in a loop
        for triangle in bottom_triangles:
            if triangle[1] in start_points:
                if top_triangles[triangle[2] - start_points[1],2] in end_points:
                    top_triangles[triangle[2] - start_points[1],2] = triangle[1]
        # Now take care of the -1's
        for start_point,end_point in zip(start_points[1:],end_points[1:]):
            if (top_triangles[start_point-start_points[1],2] == -1):
                previous_point = end_point
                while (top_triangles[previous_point-start_points[1],2] == -1):
                    previous_point = previous_point - 1
                top_triangles[start_point-start_points[1],2] = top_triangles[previous_point-start_points[1],2]
            this_point = start_point
            for this_point in range(start_point,end_point+1):
                if top_triangles[this_point-start_points[1],2] != -1:
                    last_bottom = top_triangles[this_point-start_points[1],2]
                else:
                    top_triangles[this_point-start_points[1],2] = last_bottom
        return top_triangles
        
    # First, we want a function that takes (bottom ring, top ring, top ring global index)
    #   and returns an array as long as the bottom ring array that is the closest point on the top for each bottom point
    def get_closest_points(self,bottom_ring,top_ring,global_top_ring_start_index):
        nbrs = NearestNeighbors(n_neighbors=1).fit(top_ring)
        distances, indices = nbrs.kneighbors(bottom_ring)
        return indices[:,0]  + global_top_ring_start_index
