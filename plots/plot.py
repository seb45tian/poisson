import numpy as np
from matplotlib import use
use("TkAgg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm #(best so for winter?)
import colormaps as cm #(choose from plasma, inferno, magma or viridis)

# Specify data file names
potential = "testdata.txt" #sys.argv[1]
field 	  = "field.txt"

# Read data from files
phi = np.loadtxt(infile)
E = np.loadtxt(field)

# Print out the dimensions
print(phi.shape)
N = phi.shape[0];

# Split E into its components Ex and Ey
Ex = E[0:N,:]
Ey = E[N:,:]

# Create a meshgrid
x = np.linspace(0,1.0,N)
y = np.linspace(0,1.0,N)
X, Y = np.meshgrid(x,y) # create a meshgrid from y and z coordinates!

################
##--- PLOT ---##
################

##--- UNCOMMENT FOR 3D SURFACE PLOT! ----##
# fig = plt.figure(figsize=(10,5))
# ax = fig.add_subplot(1,1,1, projection="3d")
# ax.set_title('Potential Surface' , fontsize=18, fontweight='bold')

# surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.viridis)
# ax.set_xlabel("X", fontsize=18)
# ax.set_ylabel("Y", fontsize=18)
# ax.set_zlabel(r"$\varphi$", fontsize=18)

# cb = fig.colorbar(surf, shrink=0.5, aspect=5)
# cb.set_label(label=r"$\varphi$",size=18)
# plt.tight_layout()
# plt.show()



##--- 2D FIELD PLOT! ----##

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1)
ax.set_title('Field E with contourplot of potential' , fontsize=18, fontweight='bold')
surf = ax.contourf(X, Y, phi, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.viridis)
field = ax.quiver(X, Y, Ex, Ey, pivot="middle")

ax.set_xlabel("x", fontsize=18)
ax.set_ylabel("y", fontsize=18)
cb = fig.colorbar(surf, shrink=0.5, aspect=5)
cb.set_label(label=r"$\varphi$",size=18)
plt.tight_layout()

plt.show()

