import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Problem a
alpha = [1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35]
iterations = [98, 81, 66, 53, 41, 32, 25, 45]
fig = plt.figure()
fig.set_figwidth(9)
fig.set_figheight(4)
plt.plot(alpha, iterations, '.-k', linewidth=1)
plt.xlabel(r'$\alpha_t$')
plt.ylabel(r'Number of iterations')
plt.tight_layout()
plt.grid()
fig.savefig('iterations.pdf')

# Problem b
CV = [225, 441, 625, 961, 1681]
temp = [68.195676760, 68.199187280, 68.200261791, 68.201161473, 68.201877782]
fig = plt.figure()
fig.set_figwidth(9)
fig.set_figheight(4)
plt.plot(CV, temp, '.-k', linewidth=1)
plt.xlabel(r'Number of CVs')
plt.ylabel(r'Center temperature ($^\circ$C)')
plt.ticklabel_format(useOffset=False)
plt.tight_layout()
plt.grid()
fig.savefig('center.pdf')

# Problem c
T = np.loadtxt('../cpp/solution.csv', delimiter=',')
fig = plt.figure()
fig.set_figwidth(5.5)
fig.set_figheight(4)
x, y = np.meshgrid(np.linspace(0, 0.5, num=41), np.linspace(0, 0.5, num=41))
plot = plt.pcolormesh(x, y, T)
cbar = fig.colorbar(plot)
cbar.ax.set_ylabel(r'$T$ ($^\circ$C)')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.tight_layout()
fig.savefig('solution.pdf')
