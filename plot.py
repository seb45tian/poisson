import numpy as np
from matplotlib import use
use("Qt4Agg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

infile = "testdata.txt" #sys.argv[1]

data_in = np.loadtxt(infile)

phi = data_in

N = 50;

x = np.linspace(0,1.0,N+1)
y = np.linspace(0,1.0,N+1)

u = np.zeros((N+1,N+1))

for i in range(0,N+1):
	for j in range(0,N+1):
		u[i,j] = phi[i*N+j];

################
##--- PLOT ---##
################

Y, Z = np.meshgrid(x,y) # create a meshgrid from y and z coordinates!

fig = plt.figure(figsize=(12,10))

##--- UNCOMMENT FOR 3D SURFACE PLOT! ----##
# fig.suptitle('Non-dimensional velocity profile - 2D' , fontsize=18, fontweight='bold')
ax = fig.add_subplot(1,1,1, projection="3d")

surf = ax.plot_surface(Y, Z, u, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.winter)
# ax.set_xlabel("$y^*$ - Height", fontsize=18)
# ax.set_ylabel("$Z^*$ - Width", fontsize=18)
# ax.set_zlabel("$x^*$ - Direction of flow", fontsize=18)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

##--- UNCOMMENT FOR 2D CONTOUR PLOT! ----##
# fig.suptitle('Non-dimensional velocity profile - Contour plot' , fontsize=18, fontweight='bold')
# ax = fig.add_subplot(1,1,1)
# surf = ax.contourf(Y, Z, u, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.winter)
# ax.set_xlabel("$y^*$ - Height", fontsize=18)
# ax.set_ylabel("$Z^*$ - Width", fontsize=18)
# fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

