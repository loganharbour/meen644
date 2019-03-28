import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Lx = 2
Ly = 0.02

# Load problem 4 results
u = np.loadtxt('results/u.csv', delimiter=',')
v = np.loadtxt('results/v.csv', delimiter=',')
T = np.loadtxt('results/T.csv', delimiter=',')

# Plot problem 4a
x_center = np.linspace(0, Lx, num = u.shape[1])
u_center = np.copy(u)[int(u.shape[0] / 2), :]
u_center /= u_center[0]
fig, ax = plt.subplots(1)
fig.set_figwidth(6)
fig.set_figheight(3)
ax.plot(x_center, u_center, 'k')
ax.set_xlabel(r'$x$ (m)')
ax.set_ylabel('Centerline $u / u_{{\\mathrm{{in}}}}$')
ax.grid()
fig.tight_layout()
fig.savefig('results/u_centerline.pdf', bbox_inches='tight')

# Plot problem 4b
for var, name in [(u, 'u'), (v, 'v'), (T, 'T')]:
    if name == 'T':
        dy = Ly / (var.shape[0] - 2)
        y_8 = np.hstack((0, (np.linspace(dy / 2, Ly - dy / 2, num = var.shape[0] - 2)), Ly))
    else:
        y_8 = np.linspace(0, Ly, num = var.shape[0])
    var_8 = var[:, int(0.8 * var.shape[1] / Lx)]

    fig, ax = plt.subplots(1)
    fig.set_figwidth(6)
    fig.set_figheight(3)
    if name == 'y':
        ax.semilogy(y_8, var_8, 'k')
    else:
        ax.plot(y_8, var_8, 'k')
    ax.set_xlabel(r'$y$ (m)')
    ax.set_ylabel('${}(0.8$ m$, y)$'.format(name))
    ax.grid()
    fig.tight_layout()
    fig.savefig('results/{}_0p8.pdf'.format(name), bbox_inches='tight')
