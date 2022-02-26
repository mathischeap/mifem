import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl

# create some fake data
x = y = np.arange(-4.0, 4.0, 0.02)
# here are the x,y and respective z values
X, Y = np.meshgrid(x, y)
Z = np.sinc(np.sqrt(X*X+Y*Y))
# this is the value to use for the color
V = (np.sin(Y)+1)/2

cmap = cm.bwr
mappable = cm.ScalarMappable(cmap=cmap)
mappable.set_array(V)

# create the figure, add a 3d axis, set the viewing angle
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45,60)

# here we create the surface plot, but pass V through a colormap
# to create a different color for each patch
ax.plot_surface(X, Y, Z, facecolors=cmap(V))


cbar = plt.colorbar(mappable, ax=ax, ticks=[0, 0.5, 1], #shrink=1, aspect=20,
             extend='both',
             orientation='vertical', label='Some Units')

cbar.ax.set_yticklabels(["-1", "0", "1"])  # vertically oriented colorbar

plt.show()