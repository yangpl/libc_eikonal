from discretize import TensorMesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# size of the 3D salt model:
# n1=676 n2=676 n3=210
# o1=0.0 o2=0.0 o3=0.0
# d1=20.0 d2=20.0 d3=20.0

nx = 201
ny = 201
nz = 101
dx = 10
dy = 10
dz = 5
x0 = 0
y0 = 0
z0 = -nz*dz

hx = nx*[dx,] #dx*np.ones(nx)
hy = ny*[dy,] #dy*np.ones(ny)
hz = nz*[dz,] #dz*np.ones(nz)

fvp = np.fromfile('traveltime.bin', dtype=np.float32)
vp = fvp.reshape([nz, nx, ny], order='F')
vp = np.transpose(vp, (1, 2, 0))
vp = vp[:,:,::-1] #reverse the order of z-axis


mesh = TensorMesh([hx, hy, hz], origin=(x0, y0, z0))
#print(mesh)
mesh.plot_3d_slicer(vp,  xslice=nx*dx/2, yslice=ny*dy/2, zslice=-nz*dz/2,
                xlim=(x0, (nx-1)*dx), ylim=(y0, (ny-1)*dy), zlim=(z0, 0),
                grid=[4,4,3], pcolor_opts={'cmap':'jet'})#'rainbow','jet','grey','seismic'
plt.show()



