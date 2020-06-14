from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

my_data = np.genfromtxt('cmake-build-debug/out.txt', delimiter=' ')
correct = np.genfromtxt('cmake-build-debug/correct.txt', delimiter=' ')

print("max", np.max(np.abs(my_data - correct)), np.unravel_index(np.argmax(np.abs(my_data - correct)), my_data.shape))
print("min", np.min(np.abs(my_data - correct)), np.unravel_index(np.argmin(np.abs(my_data - correct)), my_data.shape))
print("sum", np.sum(np.abs(my_data - correct)))
