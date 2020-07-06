import matplotlib.pyplot as plt 
import numpy as np

en = np.genfromtxt('energy', usecols = [1])
temp = np.genfromtxt('energy', usecols = [2])
x = np.genfromtxt('energy',usecols = [0])

plt.plot(x,en)
plt.title('Energy plot')
plt.xlabel('Timesteps dt')
plt.ylabel('Energy')
plt.show()
plt.plot(x,temp)
plt.title('Temperature plot')
plt.xlabel('Timesteps dt')
plt.ylabel('Temperature')
plt.show()
