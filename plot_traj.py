import numpy as np
import matplotlib.pyplot as plt

x = np.genfromtxt('trajectory', usecols = [1])
y = np.genfromtxt('trajectory', usecols = [2])

#plotting only the oxygen particle's trajectory
ox = x[::3]
oy = y[::3]

for i in range(len(ox)-1):
    plt.scatter(ox[i],oy[i],c='k',marker = 'o')
    plt.pause(0.001)

plt.show()