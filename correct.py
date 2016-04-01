import numpy as np
from matplotlib import use
use("Agg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# Grid
N = 51; # Use an odd number to avoid singularity?
# Domain length
xlenght = 1.0
ylength = 1.0
# create x and y vectors from 0 to x/ylength
x = np.linspace(0,xlenght,N+1)
y = np.linspace(0,ylength,N+1)
# solution matrix with 0 for phi and E
phi = np.zeros((N+1,N+1))
Ex  = []
Ey	= []
# plus and minus location
rplus  = np.array([0.3,0.3])
rminus = np.array([0.7,0.7])


# POTENTIAL: Analytical solution, breaks down at rplus and rminus! (zero division)
def potential(r):
	ret = (1/(2*np.pi))*np.log( (np.linalg.norm(r-rplus)) / (np.linalg.norm(r-rminus)) )
	return(ret)

# POTENTIAL: Analytical solution, breaks down at rplus and rminus! (zero division)
def field(r):
	ret = (1/(2*np.pi))*( (r-rplus)/((np.linalg.norm(r-rplus))**2) - (r-rminus)/((np.linalg.norm(r-rminus))**2) )
	return(ret)


# Iterate over x and y and fill solution matrix with results
for i, xi in enumerate(x):
	for j, yj in enumerate(y):
		phi[i,j] = potential(np.array([xi,yj]))
		Ex.append(field(np.array([xi,yj]))[0])
		Ey.append(field(np.array([xi,yj]))[1])



################
##--- PLOT ---##
################

X, Y = np.meshgrid(x,y) # create a meshgrid from y and z coordinates!



##--- UNCOMMENT FOR 3D SURFACE PLOT! ----##
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1, projection="3d")
ax.set_title('Potential Surface' , fontsize=18, fontweight='bold')

surf = ax.plot_surface(X, Y, phi, rstride=1, cstride=1, linewidth=0, antialiased=True, cmap = cm.viridis)
ax.set_xlabel("X", fontsize=18)
ax.set_ylabel("Y", fontsize=18)
ax.set_zlabel(r"$\varphi$", fontsize=18)

cb = fig.colorbar(surf, shrink=0.5, aspect=5)
cb.set_label(label=r"$\varphi$",size=18)
plt.tight_layout()
plt.savefig("potetential_surfaceplot_analytical.png")
# plt.show()



##--- 2D FIELD PLOT! ----##

# Normalise vector if wanted
# Ex = np.array(Ex)
# Ey = np.array(Ey)
# Ex = Ex/np.sqrt(Ex**2+Ey**2)
# Ey = Ey/np.sqrt(Ex**2+Ey**2)


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
plt.savefig("field_contourplot_analytical.png")
# plt.show()

