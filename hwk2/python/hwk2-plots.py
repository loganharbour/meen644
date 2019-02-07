import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load results
T = {}
Tinf = 25
T0 = 100
nodes = [81, 41, 21, 11, 6]
for N in nodes:
    # Remember to add Tinf to all theta and append T0 to the first
    T[N] = np.loadtxt('../results/theta_{}.csv'.format(N)) + Tinf
    T[N] = np.insert(T[N], 0, T0)

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
