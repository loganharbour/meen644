import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Lx = 2
Ly = 0.02
cp = 4183
rho = 998.3
q = 500
k = 0.609

###############################################################################
# Problem 3 prints

# Load coarse results
u = np.loadtxt('results/coarse_u.csv', delimiter=',').T
p = np.loadtxt('results/coarse_p.csv', delimiter=',').T
T = np.loadtxt('results/coarse_T.csv', delimiter=',').T
print(p.shape)

# Normalize P
p -= p[p.shape[0] - 1, :]

for var in [u, p, T]:
    for i in range(var.shape[0]):
        line = '{} & '.format(i + 1)
        for j in range(var.shape[1]):
            if var is T:
                line += '{:.5f}'.format(var[i,j])
            else:
                line += '{:.5e}'.format(var[i,j])
            if j == var.shape[1] - 1:
                line += " \\\\"
            else:
                line += " & "
        print(line)
    print()

###############################################################################
# Load refined results

u = np.loadtxt('results/fine_u.csv', delimiter=',').T
v = np.loadtxt('results/fine_v.csv', delimiter=',').T
T = np.loadtxt('results/fine_T.csv', delimiter=',').T

dx = Lx / (T.shape[0] - 2)
dy = Ly / (T.shape[1] - 2)

###############################################################################
# Plot problem 4a

x_center = np.linspace(0, Lx, num = u.shape[0])
u_center = np.copy(u)[:, int(u.shape[1] / 2)]
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

###############################################################################
# Plot problem 4b

for var, name in [(u, 'u'), (v, 'v'), (T, 'T')]:
    if name == 'T':
        y_8 = np.hstack((0, (np.linspace(dy / 2, Ly - dy / 2, num = var.shape[1] - 2)), Ly))
    else:
        y_8 = np.linspace(0, Ly, num = var.shape[1])
    var_8 = var[int(0.8 * var.shape[0] / Lx), :]

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

###############################################################################
# Plot problem 4b

# Sample u at temperature nodes
uT = np.zeros(T.shape)
for i in range(1, T.shape[0] - 1):
    for j in range(u.shape[1]):
        uT[0, j] = u[0, j]
        uT[u.shape[0], j] = u[u.shape[0] - 1, j]
        uT[i, j] = (u[i - 1, j] + u[i, j]) / 2

# Compute Nusselt numbers
Nu = np.zeros(T.shape[0])
u_left = u[0, 0]
for i in range(T.shape[0]):
    Tw = T[i, 0]
    Tb = 0
    for j in range(1, T.shape[1] - 1):
        Tb += uT[i, j] * T[i, j]
    Tb *= rho * cp * dy / (u_left * cp * Ly * rho)
    Nu[i] = 2 * Ly * q / (k * (Tw - Tb))

# Plot
x = np.hstack((0, (np.linspace(dx / 2, Lx - dx / 2, num = len(Nu) - 2)), Lx))
fig, ax = plt.subplots(1)
fig.set_figwidth(6)
fig.set_figheight(3)
ax.plot(x, Nu, 'k')
ax.set_xlabel(r'$x$ (m)')
ax.set_ylabel('Nu')
ax.grid()
fig.tight_layout()
fig.savefig('results/Nu.pdf', bbox_inches='tight')
