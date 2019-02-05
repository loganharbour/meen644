import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Prettier plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['axes.grid'] = True

def solveODEs(N, f=lambda t, y: np.array([y[1], -9.81 * np.sin(y[0])]), tf=2,
              y0=[np.pi / 9, 0], method='RK', tol=None, lamb=3):
    """Solves a set of coupled ODEs.

    Args:
        N (int): The number of timesteps.
        f (callable f(x)): The function that describes the ODEs.
        tf (float): The ending time.
        y0 (ndarray): The initial condition at t = 0.
        method (string): The solving method (RK, explicit, implicit).
        tol (float, optional): Set to enable uniform timestep refinement to
                               reach the specified tolerance.

    Returns:
        t (ndarray): The time steps the solution is solved at.
        y (ndarray): The solution vector.
    """
    # Initalize timestep size and solution storage
    dt = tf / N
    y = np.zeros((N + 1, len(y0)))
    y[0] = y0
    t = np.linspace(0, dt * N, num=N + 1)

    # Integrate over each timestep
    for n in range(N):
        if method == 'RK':
            k1 = f(t[n], y[n])
            k2 = f(t[n + 1], y[n] + dt * k1)
            y[n + 1] = y[n] + dt * (k1 + k2) / 2
        elif method == 'explicit':
            y[n + 1] = y[n] + dt * f(t[n], y[n])
        elif method == 'implicit':
            func = lambda x: y[n] + dt * f(t[n + 1], x) - x
            y[n + 1] = fsolve(func, y[n])

    # Tolerance given, need to converge
    if tol is not None:
        print('Refining in time...')
        dts, max_changes, last_changes = [], [], []
        y_old = np.copy(y)
        # Refine in time and solve again until converged
        while True:
            N *= 2
            t, y = solveODEs(N, method=method)

            # Relative change at each one-step-coarser point for each variable
            change = np.abs((y_old[1:] - y[::2][1:]) / (y_old[1:]))
            print('  dt {:.2e} s: max change = {:.2e}'.format(tf / N, np.max(change)))

            # Append for plotting
            dts.append(tf / N)
            max_changes.append(np.max(change, axis=0))
            last_changes.append(change[-1])

            # Max is less than tolerance: done
            if np.max(change) < tol:
                return t, y, N, dts, np.array(max_changes), np.array(last_changes)

            # Didn't converge, copy previous solution
            y_old = np.copy(y)

    # Standard run without convergence, return just t and y
    return t, y

################################################################################
# Part a: RK solution with dt = 0.5 (N = 4) and plot

t_coarse, y_coarse = solveODEs(4)

fig, ax = plt.subplots(1, 2)
fig.set_figwidth(9)
fig.set_figheight(3)
plt.rcParams['lines.linewidth'] = 1.0
ax[0].plot(t_coarse, y_coarse[:,0], '--k.')
ax[0].set_ylabel(r'$\theta$')
ax[0].set_xlabel(r'$t$ (sec)')
ax[1].plot(t_coarse, y_coarse[:,1], '--k.')
ax[1].set_ylabel(r'$\omega$')
ax[1].set_xlabel(r'$t$ (sec)')
fig.savefig('a.pdf', bbox_inches='tight')

################################################################################
# Part b: Determine grid independence and plot

t_rk, y_rk, N_converged, dts, max_changes, last_changes = solveODEs(4, tol=1e-4)
dt_converged = dts[-1]

fig = plt.figure()
fig.set_figwidth(9)
fig.set_figheight(4)
plt.loglog(dts, max_changes[:,0], 'k.--', label=r'$\theta$: all coarse timesteps')
plt.loglog(dts, last_changes[:,0], 'k.:', label=r'$\theta: t = 2$ s')
plt.loglog(dts, max_changes[:,1], 'b.--', label=r'$\omega$: all coarse timesteps')
plt.loglog(dts, last_changes[:,1], 'b.:', label=r'$\omega: t = 2$ s')
plt.xlabel(r'$\Delta t$ (s)')
plt.ylabel('Relative change')
plt.tight_layout()
plt.legend()
fig.savefig('b.pdf')

################################################################################
# Part c: All three methods

t_exp, y_exp = solveODEs(N_converged, method='explicit')
t_imp, y_imp = solveODEs(N_converged, method='implicit')

fig, ax = plt.subplots(1, 2)
fig.set_figwidth(9)
fig.set_figheight(3)
plt.rcParams['lines.linewidth'] = 2.0
ax[0].plot(t_exp, y_exp[:,0], label='Explicit RK')
ax[0].plot(t_exp, y_exp[:,0], '--', label='Explicit Euler')
ax[0].plot(t_imp, y_imp[:,0], ':', label='Implicit Euler')
ax[0].set_ylabel(r'$\theta$')
ax[0].set_xlabel(r'$t$ (sec)')
ax[1].plot(t_exp, y_exp[:,1])
ax[1].plot(t_exp, y_exp[:,1], '--')
ax[1].plot(t_imp, y_imp[:,1], ':')
ax[1].set_ylabel(r'$\omega$')
ax[1].set_xlabel(r'$t$ (sec)')
handles, labels = ax[0].get_legend_handles_labels()
lgd = ax[1].legend(handles, labels, loc='lower center', bbox_to_anchor=(-0.2, -0.34), ncol=3, fontsize=9)
fig.savefig('c.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
