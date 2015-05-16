import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

header = np.loadtxt('header')
ndim = header[0]
saved_rings = header[1]
steps_per_save = header[2]
dt = header[3]

rings = np.loadtxt('output')

fig = plt.figure()
ax = fig.gca(projection='3d')

for i in range(24,26):
    myrows = np.where( rings[:,0] == i )
    ax.plot(rings[myrows,1][0],rings[myrows,2][0],rings[myrows,3][0])
#ax.plot(rings[-1,:,0],rings[-1,:,1],rings[-1,:,2])
plt.show()
