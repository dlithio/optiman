import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

dims = [1,2,4]

header = np.loadtxt('header')
ndim = header[0]
saved_rings = header[1]
steps_per_save = header[2]
dt = header[3]

def read_array(filename,datatype):
    try:
        result = np.fromfile(filename,dtype=datatype)
    except IOError:
        print(filename + " was not yet available")
    return result
        
rings = read_array('output',np.float64)

numrows = rings.shape[0]/(ndim+1)
rings.shape = (numrows,ndim+1)

fig = plt.figure()
ax = fig.gca(projection='3d')

for i in range(1,int(saved_rings+1)):
    myrows = np.where( rings[:,0] == i )
    ax.plot(rings[myrows,dims[0]][0],rings[myrows,dims[1]][0],rings[myrows,dims[2]][0])
#ax.plot(rings[-1,:,0],rings[-1,:,1],rings[-1,:,2])
plt.show()
