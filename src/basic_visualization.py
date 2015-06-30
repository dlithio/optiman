import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# Control which steady states are plotted
ss = ['9','23','37','59','71','85','97','110','123','133','143']


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
       
##Load files
#rings = read_array('output_59',np.float64)
#numrows = rings.shape[0]/(ndim+1)
#rings.shape = (numrows,ndim+1)

rings2 = read_array('output',np.float64)
fdot = read_array('fdot',np.float64)
numrows = rings2.shape[0]/(ndim+1)
rings2.shape = (numrows,ndim+1)

# Plot files
# Controls which dimensions are plotted
dims = [1,2,4]
#
fig = plt.figure()
ax = fig.gca(projection='3d')

#saved_rings = 25 

#for i in range(1,int(saved_rings+1)):
    #myrows = np.where( rings[:,0] == i )
    #ax.plot(rings[myrows,dims[0]][0],rings[myrows,dims[1]][0],rings[myrows,dims[2]][0])
##ax.plot(rings[-1,:,0],rings[-1,:,1],rings[-1,:,2])

for i in range(1,int(saved_rings+1)):
    myrows = np.where( rings2[:,0] == i )
    ax.plot(rings2[myrows,dims[0]][0],rings2[myrows,dims[1]][0],rings2[myrows,dims[2]][0])
#ax.plot(rings[-1,:,0],rings[-1,:,1],rings[-1,:,2])

for pointname in ss:
    newpoint = np.loadtxt(pointname)
    ax.scatter(newpoint[dims[0]-1],newpoint[dims[1]-1],newpoint[dims[2]-1])

kx = np.loadtxt("kx_projections")
ky = np.loadtxt("ky_projections")
ax.set_xlabel("Mode "+str(dims[0])+"=("+str(int(kx[dims[0]-1]))+","+str(int(ky[dims[0]-1]))+")")
ax.set_ylabel("Mode "+str(dims[1])+"=("+str(int(kx[dims[1]-1]))+","+str(int(ky[dims[1]-1]))+")")
ax.set_zlabel("Mode "+str(dims[2])+"=("+str(int(kx[dims[2]-1]))+","+str(int(ky[dims[2]-1]))+")")
plt.show()


#i = 26
#myrows = np.where( rings2[:,0] == i )
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(rings2[myrows,dims[0]][0][::4],rings2[myrows,dims[1]][0][::4],rings2[myrows,dims[2]][0][::4],c=((fdot[myrows,1][0][::4]+1.0)/2.0),cmap='bwr')
#plt.show()
