from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

my_data = np.genfromtxt('cmake-build-debug/out.txt', delimiter=' ')
correct = np.genfromtxt('cmake-build-debug/correct.txt', delimiter=' ')

print(np.sum(my_data - correct))
