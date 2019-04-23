import numpy as np
import matplotlib.pyplot as plt
import glob
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Lx = 0.6
Ly = 0.02
b = 0.06
q = 16
rho = 998.3
k = 0.609
cp = 4183
mu = 0.001002

###############################################################################
# Problem 1

Ny_vals = [30, 50, 70, 90]
Nx = 160
dx = Lx / Nx

# Reattachment lengths, u sampled at T points
xr, T, Nu_top, Nu_bot = {}, {}, {}, {}
for Ny in Ny_vals:
    dy = Ly / Ny
    u = np.loadtxt('results/Re200_Nx160_Ny{}_u.csv'.format(Ny), delimiter=',').T
    T[Ny] = np.loadtxt('results/Re200_Nx160_Ny{}_T.csv'.format(Ny), delimiter=',').T
    # Reattachment lengths
    u_top = u[:, Ny]
    for i in range(int(np.floor(b / dx)) + 1, Nx):
        if u_top[i] > 0:
            m = (u_top[i] - u_top[i - 1]) / dx
            xr[Ny] = i * dx - (u_top[i] / m) - b
            break
    # u at T points
    uT = np.zeros(T[Ny].shape)
    for i in range(1, T[Ny].shape[0] - 1):
        for j in range(u.shape[1]):
            uT[0, j] = u[0, j]
            uT[u.shape[0], j] = u[u.shape[0] - 1, j]
            uT[i, j] = (u[i - 1, j] + u[i, j]) / 2
    # Nusselt numbers
    Nu_bot[Ny] = []
    Nu_top[Ny] = []
    T_bulk_denom = 0
    for j in range(u.shape[1]):
        T_bulk_denom += cp * dy * u[0, j] * rho
        # T_bulk_denom += cp * dy * uleft * rho
    for i in range(T[Ny].shape[0]):
        T_bulk = 0
        for j in range(1, T[Ny].shape[1] - 1):
            if i * dx <= b and j * dy >= Ly / 2:
                continue
            T_bulk += rho * uT[i, j] * cp * T[Ny][i, j] * dy
        T_bulk /= T_bulk_denom
        T_bot, T_top = T[Ny][i, 0], T[Ny][i, Ny + 1]
        Nu_bot[Ny].append(Ly * q / (k * (T_bot - T_bulk)))
        if i > b / dx:
            Nu_top[Ny].append(Ly * q / k * (T_top - T_bulk))

print('Problem 2 reattachment (160xNy grid, Re = 200):')
print('  (Ny: xr):', xr)

# Problem 1a
dy = Ly / 90
T_y = np.hstack((0, (np.linspace(dy / 2, Ly - dy / 2, num = T[90].shape[1] - 2)), Ly))
fig, ax = plt.subplots(1)
fig.set_figwidth(8)
fig.set_figheight(3)
for x in (Ly * np.array([6, 12, 24]) + b):
    i_mid = x / dx + 0.5
    i_min, i_max = int(np.floor(i_mid)), int(np.ceil(i_mid))
    ax.plot(T_y, (T[90][i_min, :] + T[90][i_max, :]) / 2, label='x = {:.2f} m'.format(x), linewidth=1)
plt.xlabel('$y$ (m)')
plt.ylabel('$T(x, y)$')
plt.legend()
plt.grid()
plt.tight_layout()
fig.savefig('results/1a_Ty.pdf')

# Problem 1b
T_x = np.hstack((0, (np.linspace(dx / 2, Lx - dx / 2, num = T[90].shape[0] - 2)), Lx))
fig, ax = plt.subplots(1)
fig.set_figwidth(8)
fig.set_figheight(3)
after_step = np.where(T_x >= b)[0][0]
ax.plot(T_x, T[90][:, 0], label='Bottom wall', linewidth=1)
ax.plot(T_x[after_step:], T[90][after_step:, -1], label='Top wall'.format(Ly), linewidth=1)
plt.xlabel('$x$ (m)')
plt.ylabel('$T(x)$')
plt.legend()
plt.grid()
plt.tight_layout()
fig.savefig('results/1b_Twall.pdf')

# Problem 1c
fig, ax = plt.subplots(1, 2)
fig.set_figwidth(8)
fig.set_figheight(3)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
i = 0
for Ny in Ny_vals:
    ax[0].plot(T_x[np.where(T_x > b)], Nu_top[Ny], label='$N_y = {}$'.format(Ny), linewidth=1)
    ax[1].plot(T_x, Nu_bot[Ny], linewidth=1)
ax[0].set_xlabel('$x$ (m)')
ax[1].set_xlabel('$x$ (m)')
ax[0].set_ylabel('Nu$(x)$')
ax[0].set_title('Top wall', fontsize=10)
ax[1].set_title('Bottom wall', fontsize=10)
ax[0].grid()
ax[1].grid()
handles, labels = ax[0].get_legend_handles_labels()
lgd = ax[0].legend(handles, labels, loc='lower center', bbox_to_anchor=(1.0, -0.37),
                ncol=4, fontsize=9)
fig.tight_layout()
fig.savefig('results/1c_Nu.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))

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
    # Cut along various xcome
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
fig.set_figheight(3)
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
lgd = ax[1].legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.37),
                ncol=3, fontsize=9)
fig.tight_layout()
fig.savefig('results/2_u.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))

# u velocity profile
fig, ax = plt.subplots(1, 3)
fig.set_figwidth(8)
fig.set_figheight(3)
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
lgd = ax[1].legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.37),
                ncol=3, fontsize=9)
fig.tight_layout()
fig.savefig('results/2_v.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
