import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load results
T = {}
nodes = [6, 11, 21, 41, 81]
for N in nodes:
    T[N] = np.loadtxt('../results/result_{}.csv'.format(N))

# Plot results
fig = plt.figure()
fig.set_figwidth(9)
fig.set_figheight(4)
for N in nodes:
    x = np.linspace(0, 1, num=N)
    plt.plot(x, T[N], '.-', linewidth=1, label='ITMAX = {}'.format(N))
plt.xlabel(r'$x$ (m)')
plt.ylabel(r'$T$ ($^\circ$C)')
plt.tight_layout()
plt.grid()
plt.legend()
fig.savefig('result.pdf')
