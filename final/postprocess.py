import numpy as np
import matplotlib.pyplot as plt
import glob
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Lx = 0.6
Ly = 0.02
b = 0.06

###############################################################################
# Problem 1

# Reattachment lengths for 160xNy
Ny_vals = [30, 50, 70, 90]
Nx = 160
dx = Lx / Nx
xr = {}
for Ny in Ny_vals:
    u_top = np.loadtxt('results/Re200_Nx160_Ny{}_u.csv'.format(Ny), delimiter=',').T[:, Ny]
    for i in range(int(np.floor(b / dx)) + 1, Nx):
        if u_top[i] > 0:
            m = (u_top[i] - u_top[i - 1]) / dx
            xr[Ny] = i * dx - (u_top[i] / m) - b
            break
print('Problem 2 reattachment (160xNy grid, Re = 200):')
print('  (Ny: xr):', xr)

# Data for 1a-1b
T_data = np.loadtxt('results/Re200_Nx160_Ny90_T.csv', delimiter=',').T
Ny = 90
dy = Ly / Ny
T_x = np.hstack((0, (np.linspace(dx / 2, Lx - dx / 2, num = T_data.shape[0] - 2)), Lx))
T_y = np.hstack((0, (np.linspace(dy / 2, Ly - dy / 2, num = T_data.shape[1] - 2)), Ly))

# Problem 1a
fig, ax = plt.subplots(1)
fig.set_figwidth(6)
fig.set_figheight(3)
for x in (Ly * np.array([6, 12, 24]) + b):
    i_mid = x / dx + 0.5
    i_min, i_max = int(np.floor(i_mid)), int(np.ceil(i_mid))
    T = (T_data[i_min, :] + T_data[i_max, :]) / 2
    ax.plot(T_y, T, label='x = {:.2f} m'.format(x), linewidth=1)
plt.xlabel('$y$ (m)')
plt.ylabel('$T(x, y)$')
plt.legend()
plt.grid()
plt.tight_layout()
fig.savefig('results/1a_Ty.pdf')

# Problem 1b
fig, ax = plt.subplots(1)
fig.set_figwidth(6)
fig.set_figheight(3)
after_step = np.where(T_x >= b)[0][0]
ax.plot(T_x, T_data[:, 0], label='y = 0 m', linewidth=1)
ax.plot(T_x[after_step:], T_data[after_step:, -1], label='y = {} m'.format(Ly), linewidth=1)
plt.xlabel('$x$ (m)')
plt.ylabel('$T(x, y)$')
plt.legend()
plt.grid()
plt.tight_layout()
fig.savefig('results/1b_Twall.pdf')

###############################################################################
# Problem 2

Ny = 70
dx = Lx / 160
dy = Ly / Ny

x_vals = (Ly * np.array([6, 12, 24]) + b)
Re_vals = [100, 300, 400]
xr = {}
u_cut, v_cut = {}, {}
for Re in Re_vals:
    u = np.loadtxt('results/Re{}_Nx160_Ny70_u.csv'.format(Re), delimiter=',').T
    v = np.loadtxt('results/Re{}_Nx160_Ny70_v.csv'.format(Re), delimiter=',').T
    # Cut along various x
    u_cut[Re], v_cut[Re] = {}, {}
    for x in x_vals:
        u_i = x / dx
        v_i = x / dx + 0.5
        u_i_min, u_i_max = int(np.ceil(u_i)), int(np.ceil(u_i))
        v_i_min, v_i_max = int(np.ceil(v_i)), int(np.ceil(v_i))
        u_cut[Re][x] = (u[u_i_min, :] + u[u_i_min, :]) / 2
        v_cut[Re][x] = (v[u_i_min, :] + v[u_i_min, :]) / 2
    # Reattachment length
    u_top = u[:, Ny]
    for i in range(int(np.floor(b / dx)) + 1, Nx):
        if u_top[i] > 0:
            m = (u_top[i] - u_top[i - 1]) / dx
            xr[Re] = i * dx - (u_top[i] / m) - b
            break
print('Problem 2 reattachment (160x70 grid):')
print('  (Re: xr):', xr)

u_y = np.hstack((0, (np.linspace(dy / 2, Ly - dy / 2, num = len(u_cut[Re][x]) - 2)), Ly))
v_y = np.linspace(0, Ly, num=len(v_cut[Re][x]))

# u velocity profile
fig, ax = plt.subplots(1, 3)
fig.set_figwidth(8)
fig.set_figheight(2.5)
i = 0
for Re in Re_vals:
    for x in x_vals:
        ax[i].plot(u_y, u_cut[Re][x], label='x = {:.2f} m'.format(x), linewidth=1)
    ax[i].set_title('Re = {}'.format(Re), fontsize=10)
    ax[i].set_xlabel('$y$ (m)')
    ax[i].grid()
    i += 1
ax[0].set_ylabel('$u(x, y)$')
handles, labels = ax[1].get_legend_handles_labels()
lgd = ax[1].legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.48),
                ncol=3, fontsize=9)
fig.tight_layout()
fig.savefig('results/2_u.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))

# u velocity profile
fig, ax = plt.subplots(1, 3)
fig.set_figwidth(8)
fig.set_figheight(2.5)
i = 0
for Re in Re_vals:
    for x in x_vals:
        ax[i].plot(v_y, v_cut[Re][x], label='x = {:.2f} m'.format(x), linewidth=1)
    ax[i].set_title('Re = {}'.format(Re), fontsize=10)
    ax[i].set_xlabel('$y$ (m)')
    ax[i].grid()
    i += 1
ax[0].set_ylabel('$v(x, y)$')
handles, labels = ax[1].get_legend_handles_labels()
lgd = ax[1].legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.48),
                ncol=3, fontsize=9)
fig.tight_layout()
fig.savefig('results/2_v.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
