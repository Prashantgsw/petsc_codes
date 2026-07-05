import numpy as np
import matplotlib.pyplot as plt

# Ask for the data file
filename = input("Enter data file name: ")
data = np.loadtxt(filename)

x    = data[:,0]
rho  = data[:,1]
rhou = data[:,2]
rhoE = data[:,3]

plt.figure()
plt.plot(x, rho)
plt.xlabel("x")
plt.ylabel(r"$\rho$")
plt.grid(True)

plt.figure()
plt.plot(x, rhou)
plt.xlabel("x")
plt.ylabel(r"$\rho u$")
plt.grid(True)

plt.figure()
plt.plot(x, rhoE)
plt.xlabel("x")
plt.ylabel(r"$\rho E$")
plt.grid(True)

plt.show()