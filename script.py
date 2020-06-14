from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

X = np.arange(-100, 101, 1)
Y = np.arange(-100, 101, 1)


X, Y = np.meshgrid(X, Y)


my_data = np.genfromtxt('cmake-build-debug/out.txt', delimiter=',')
dif_data = np.genfromtxt('cmake-build-debug/out2.txt', delimiter=',')

plt.figure()
plt.imshow(my_data)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, my_data, cmap=cm.coolwarm,
                 linewidth=0, antialiased=False)
ax.set_zlim(-.000001, .000001)
plt.show()
