import numpy as np
import matplotlib.pyplot as plt
import glob
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Problem a
fig = plt.figure()
fig.set_figwidth(9)
fig.set_figheight(4)
alphas = [1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35]
for alpha in alphas:
    residuals = np.loadtxt('results/a/w{}.csv'.format(alpha))
    iterations = list(range(1, len(residuals) + 1))
    plt.semilogy(iterations, residuals, '-', linewidth=1,
                 label=r'$\alpha = ${:.2f}'.format(alpha))
plt.xlabel(r'Iteration count')
plt.ylabel(r'Residual')
plt.legend()
plt.tight_layout()
plt.grid()
fig.savefig('results/a.pdf')

# Problem b number of iterations
fig, ax = plt.subplots(2, 2)
fig.set_figwidth(11)
fig.set_figheight(8)
alphas = [1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
Ns = [[21, 25], [31, 41]]
for i in range(2):
    for j in range(2):
        N = Ns[i][j]
        for alpha in alphas:
            residuals = np.loadtxt('results/b-iterations/N{}_w{}.csv'.format(N, alpha))
            iterations = list(range(1, len(residuals) + 1))
            ax[i][j].semilogy(iterations, residuals, '-', linewidth=1,
                              label=r'$\alpha = ${:.2f}'.format(alpha))
        ax[i][j].set_title('{} CVs'.format(N * N))
        ax[i][j].grid()
        if i == 1: ax[i][j].set_xlabel('Iteration count')
        if j == 0: ax[i][j].set_ylabel('Residual')
handles, labels = ax[1][1].get_legend_handles_labels()
lgd = ax[1][1].legend(handles, labels, loc='lower center', bbox_to_anchor=(-0.2, -0.34),
                      ncol=4, fontsize=9)
fig.savefig('results/b-iterations.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))

# Problem b temps
fig = plt.figure()
fig.set_figwidth(9)
fig.set_figheight(4)
Ns = [15**2, 21**2, 26**2, 31**2, 41**2]
for file in glob.glob('results/b-temps/*.csv'):
    alpha = float(os.path.splitext(file)[0].split("/")[-1])
    temps = np.loadtxt(file)
    plt.plot(Ns, temps, '.-', linewidth=1, label=r'$\alpha =$ {:.2f}'.format(alpha))
plt.xlabel(r'Number of CVs')
plt.ylabel(r'Center plate temperature ($^\circ$C)')
plt.ticklabel_format(useOffset=False)
plt.legend()
plt.tight_layout()
plt.grid()
fig.savefig('results/b-temps.pdf')

# Problem c
T = np.loadtxt('results/c/solution.csv', delimiter=',')
fig = plt.figure()
fig.set_figwidth(5.5)
fig.set_figheight(4)
x, y = np.meshgrid(np.linspace(0, 0.5, num=41), np.linspace(0, 0.5, num=41))
plot = plt.contourf(x, y, T, levels=np.linspace(50, 100, num=25))
cbar = fig.colorbar(plot)
cbar.set_ticks(np.arange(50, 101, 10))
cbar.ax.set_ylabel(r'$T$ ($^\circ$C)')
plt.xlabel('$x$ (m)')
plt.ylabel('$y$ (m)')
plt.tight_layout()
fig.savefig('results/c.pdf')
