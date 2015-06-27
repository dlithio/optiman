import numpy as np
# For the scatter plot
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pdb
from sklearn.neighbors import NearestNeighbors
# For mayavi
from tvtk.api import tvtk
from mayavi.scripts import mayavi2

# Set up the example points
num_bands=4
rings = np.zeros((16+16+32+16,num_bands))
# Set the rings numbers
rings[:16,0] = 1.0
rings[16:32,0] = 2.0
rings[32:64,0] = 3.0
rings[64:,0] = 4.0
# Set the z coordinates
rings[:16,3] = 1.0
rings[16:32,3] = 2.0
rings[32:64,3] = 3.0
rings[64:,3] = 4.0
# Set the roll for each level
roll2 = -3
roll3 = 4
roll4 = 11
# Set the x coordinates
n=16.0
rings[:16,1] = np.sin(np.linspace(0,2*np.pi-2*np.pi/n,n))
rings[16:32,1] = np.roll(np.sin(np.linspace(0,2*np.pi-2*np.pi/n,n)),roll2)
rings[64:,1] = np.roll(np.sin(np.linspace(0,2*np.pi-2*np.pi/n,n)),roll4)
n = 32
rings[32:64,1] = np.roll(np.sin(np.linspace(0,2*np.pi-2*np.pi/n,n)),roll3)
# Set the y coordinates
n=16
rings[:16,2] = np.cos(np.linspace(0,2*np.pi-2*np.pi/n,n))
rings[16:32,2] = np.roll(np.cos(np.linspace(0,2*np.pi-2*np.pi/n,n)),roll2)
rings[64:,2] = np.roll(np.cos(np.linspace(0,2*np.pi-2*np.pi/n,n)),roll4)
n = 32
rings[32:64,2] = np.roll(np.cos(np.linspace(0,2*np.pi-2*np.pi/n,n)),roll3)

## Get a view to make sure the points are what we want
## Plot files
## Controls which dimensions are plotted
#dims = [1,2,3]
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#saved_rings = num_bands
#for i in range(1,int(saved_rings+1)):
    #myrows = np.where( rings[:,0] == i )
    #ax.scatter(rings[myrows,dims[0]][0],rings[myrows,dims[1]][0],rings[myrows,dims[2]][0])
#ax.set_xlabel("x")
#ax.set_ylabel("y")
#ax.set_zlabel("z")
#plt.show()

#tester_manifold = Manifold.Manifold(rings.shape[0],num_bands)
# Now the test of the new visualization software
def read_array(filename,datatype):
    try:
        result = np.fromfile(filename,dtype=datatype)
    except IOError:
        print(filename + " was not yet available")
    return result
    
# Load header for loading real file
header = np.loadtxt('header')
ndim = header[0]
saved_rings = header[1]
steps_per_save = header[2]
dt = header[3]
# Load an actual rings file
rings2 = read_array('output',np.float64)
numrows = rings2.shape[0]/(ndim+1)
rings2.shape = (numrows,ndim+1)
dims = [1,2,4]
full_dims = [0]
full_dims.extend(dims)
#myrows = np.where( (rings2[:,0] <= 17) and (rings2[:,0] >= 16))
rings = rings2[(rings2[:,0] >= 1) & (rings2[:,0] <= 10000)][:,full_dims]

# Function to get start points of each ring
def get_start_points(rings):
    start_points = []
    ring_num = 0
    for row_num,row_ring_num in enumerate(rings[:,0]):
        if row_ring_num != ring_num:
            ring_num = row_ring_num
            start_points.append(row_num)
    return start_points

start_points = get_start_points(rings)
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

# First, we want a function that takes (bottom ring, top ring, top ring global index)
#   and returns an array as long as the bottom ring array that is the closest point on the top for each bottom point
def get_closest_points(bottom_ring,top_ring,global_top_ring_start_index):
    nbrs = NearestNeighbors(n_neighbors=1).fit(top_ring)
    distances, indices = nbrs.kneighbors(bottom_ring)
    return indices[:,0]  + global_top_ring_start_index
# Use the function to build the bottom triangles
for ring_num in range(len(start_points) - 1):
    bottom_ring_start = start_points[ring_num]
    bottom_ring_end = end_points[ring_num] + 1
    top_ring_start = start_points[ring_num + 1]
    top_ring_end = end_points[ring_num + 1] + 1
    bottom_ring = rings[bottom_ring_start:bottom_ring_end,1:]
    top_ring = rings[top_ring_start:top_ring_end,1:]
    bottom_triangles[bottom_ring_start:bottom_ring_end,2] = \
                 get_closest_points(bottom_ring,top_ring,top_ring_start)

# Next, we want a function that takes (bottom ring, top ring, closest points to bottom on the top)
#   and returns an array as long as the top ring array that is the point on the bottom ring that we
#   should connect to when the corresponding point in the top is first in the triangle
def get_top_to_bottom(bottom_triangles,top_triangles,start_points,end_points):
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
get_top_to_bottom(bottom_triangles,top_triangles,start_points,end_points)


# Now visualize it
# The TVTK dataset.
all_triangles = np.concatenate((bottom_triangles,top_triangles), axis=0)
mesh = tvtk.PolyData(points=rings[:,1:], polys=all_triangles)
#mesh = tvtk.PolyData(points=rings[:,1:], polys=bottom_triangles[:2,:])
mesh.point_data.scalars = rings[:,0]
mesh.point_data.scalars.name = 'Time'
# The function to make mayavi work
@mayavi2.standalone
def view():
    from mayavi.sources.vtk_data_source import VTKDataSource
    from mayavi.modules.surface import Surface
    
    mayavi.new_scene()
    src = VTKDataSource(data = mesh)
    mayavi.add_source(src)
    s = Surface()
    mayavi.add_module(s)
view()

#for i in range(1,int(saved_rings+1)):
    #myrows = np.where( rings[:,0] == i )
    #temp_ring = np.zeros((myrows[0].shape[0],3))
    #temp_ring[:,0] = rings[myrows,dims[0]][0]
    #temp_ring[:,1] = rings[myrows,dims[1]][0]
    #temp_ring[:,2] = rings[myrows,dims[2]][0]
    #tester_manifold.add_ring(temp_ring)
#tester_manifold.make_triangles()
#tester_manifold.set_colors()
#my_points = tester_manifold.points[:tester_manifold.current_row,:]
#ax.scatter(xs=my_points[:,0],ys=my_points[:,1],zs=my_points[:,2])
#plt.show()
