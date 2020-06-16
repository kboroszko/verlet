from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from matplotlib import cm

def f(X, a,b):
    return (a/X) + b

my_data = np.genfromtxt('results.txt', delimiter=',')

param, param_cov = curve_fit(f, my_data[:,0], my_data[:,1])
ans = f(my_data[:,0], *param)

total = param[0] + param[1]
s = param[1]/total
p = param[0]/total

print('s', param[1])
print('p', param[0])
print('max speedup', 1/s)

N = my_data[:,0]


plt.plot(N, ans, '--', color ='red', label =("fit p={} s={}".format(param[0], param[1])))
plt.scatter(N, my_data[:,1], label='data')
plt.legend()
plt.xlabel('N proc')
plt.ylabel('T [us]')
plt.title('data fit t=s+p/N')

plt.figure()
plt.plot(N, 1./(s + p/N), color ='red' )
plt.xlabel('N proc')
plt.ylabel('speedup')
plt.title('hard speedup (Amdahl’s law)')

plt.figure()
plt.plot(N, (s + p*N), color ='red' )
plt.xlabel('N proc')
plt.ylabel('speedup')
plt.title('soft speedup (Gustafson’s law )')
plt.show()




