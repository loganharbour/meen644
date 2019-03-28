import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Lx = 2

# Load problem 4 results
u = np.loadtxt('results/u.csv', delimiter=',')
v = np.loadtxt('results/v.csv', delimiter=',')
T = np.loadtxt('results/T.csv', delimiter=',')

# Problem 4 centerline velocity
x_center = np.linspace(0, Lx, num = u.shape[1])
u_center = u[int(u.shape[0] / 2), :]
u_center /= u_center[0]
