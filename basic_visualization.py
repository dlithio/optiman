import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

with open('output','r') as f:
    first_line = f.readline()

first_line = first_line.replace('\n','')
first_line = first_line.split()
ndim = int(first_line[0])
npoints = int(first_line[1])

rings = np.loadtxt('output',skiprows=1)

rings.shape = (rings.shape[0],npoints,ndim)

fig = plt.figure()
ax = fig.gca(projection='3d')

#for i in range(rings.shape[0]):
    #ax.plot(rings[i,:,0],rings[i,:,1],rings[i,:,2])
ax.plot(rings[-1,:,0],rings[-1,:,1],rings[-1,:,2])
plt.show()
