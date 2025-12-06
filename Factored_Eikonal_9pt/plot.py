import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from discretize import TensorMesh
import numpy as np
from matplotlib.colors import LogNorm

# size of the 3D salt model:
# n1=676 n2=676 n3=210
# o1=0.0 o2=0.0 o3=0.0
# d1=20.0 d2=20.0 d3=20.0

nx = 101
ny = 101
nz = 51
dx = 100
dy = 100
dz = 100
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

# 直接美化当前图形的函数
def add_axis_labels():
    """添加坐标轴标签和单位到当前图形"""
    # 获取当前图形的所有轴
    axs = plt.gcf().get_axes()
    
    # 假设第一个轴是XY切片，第二个是XZ切片，第三个是ZY切片，第四个是颜色条
    if len(axs) >= 4:
        # XY切片
        axs[0].set_xlabel("x (m)")
        axs[0].set_ylabel("y (m)")
        
        # XZ切片
        axs[1].set_xlabel("x (m)")
        axs[1].set_ylabel("z (m)")
        
        # ZY切片
        # axs[2].set_xlabel("Depth (m)")
        # axs[2].set_ylabel("Northing (m)")
        
        # 颜色条
        axs[3].set_ylabel("Traveltime (s)")

mesh = TensorMesh([hx, hy, hz], origin=(x0, y0, z0))

# 直接调用 plot_3d_slicer，不创建新图形
mesh.plot_3d_slicer(
    vp,  
    xslice=nx*dx/2, 
    yslice=ny*dy/2, 
    zslice=-nz*dz/2,
    xlim=(x0, (nx-1)*dx), 
    ylim=(y0, (ny-1)*dy), 
    zlim=(z0, 0),
    grid=[4,4,3], 
    pcolor_opts={'cmap':'jet'}
)

# 添加坐标轴标签和单位
add_axis_labels()

# 显示图形
plt.show()
