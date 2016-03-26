import numpy as np
from matplotlib import use
use("TkAgg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm #(best so for winter?)
import colormaps as cm #(choose from plasma, inferno, magma or viridis)

infile = "testdata.txt" #sys.argv[1]
# infile = sys.argv[1]
# outfile = sys.argv[2]


data_in = np.loadtxt(infile)

phi = data_in

print(phi.shape)


N = phi.shape[0];

x = np.linspace(0,1.0,N)
y = np.linspace(0,1.0,N)

################
##--- PLOT ---##
################

X, Y = np.meshgrid(x,y) # create a meshgrid from y and z coordinates!

##--- UNCOMMENT FOR 3D SURFACE PLOT! ----##
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,1,1, projection="3d")
ax.set_title('Potential Surface' , fontsize=18, fontweight='bold')

surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.viridis)
ax.set_xlabel("X", fontsize=18)
ax.set_ylabel("Y", fontsize=18)
ax.set_zlabel(r"$\varphi$", fontsize=18)

cb = fig.colorbar(surf, shrink=0.5, aspect=5)
cb.set_label(label=r"$\varphi$",size=18)
plt.tight_layout()
plt.show()



##--- 2D FIELD PLOT! ----##

# Normalise vector if wanted
# Ex = np.array(Ex)
# Ey = np.array(Ey)
# Ex = Ex/np.sqrt(Ex**2+Ey**2)
# Ey = Ey/np.sqrt(Ex**2+Ey**2)
# fig = plt.figure(figsize=(10,8))
# ax = fig.add_subplot(1,1,1)
# ax.set_title('Field E with contourplot of potential' , fontsize=18, fontweight='bold')
# surf = ax.contourf(X, Y, phi, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.viridis)
# field = ax.quiver(X, Y, Ex, Ey, pivot="middle")

# ax.set_xlabel("x", fontsize=18)
# ax.set_ylabel("y", fontsize=18)
# cb = fig.colorbar(surf, shrink=0.5, aspect=5)
# cb.set_label(label=r"$\varphi$",size=18)
# plt.tight_layout()

# plt.show()

