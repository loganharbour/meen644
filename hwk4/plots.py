import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load problem 2 results
L = 0.2
Re = 400
rho = 998.3
mu = 1.002e-3
u_ref = Re * mu / (rho * L)
Ns = [8, 16, 32, 64, 128, 256]
x, y, u, v = {}, {}, {}, {}
for N in Ns:
    valu = np.loadtxt('results/{}_u.csv'.format(N), delimiter=',')
    u[N] = valu[:, int(N / 2)] / u_ref
    y[N] = np.linspace(0, 1, num=len(u[N]))
    valv = np.loadtxt('results/{}_v.csv'.format(N), delimiter=',')
    v[N] = valv[int(N / 2), :] / u_ref
    x[N] = np.linspace(0, 1, num=len(v[N]))

# Problem 2 reference
ref_u = [0, -0.03177, -0.06189, -0.08923, -0.11732, -0.14410, -0.17257, -0.19996,
         -0.22849, -0.25458, -0.27956, -0.11446, 0.35427, 0.37502, 0.40187,
         0.43679, 0.48695, 0.54512, 0.61626, 0.70001, 0.79419, 0.90472, 1.0]
ref_y = [0, 0.02, 0.0405, 0.0601, 0.0806, 0.1001, 0.1206, 0.1401, 0.1606, 0.1802,
         0.2007, 0.5005, 0.9009, 0.9106, 0.9204, 0.9302, 0.9409,  0.9507, 0.9604,
         0.9702, 0.9800, 0.9907, 1]
ref_v = [0, 0.05951, 0.11028, 0.14906, 0.18047, 0.20578, 0.22746, 0.24397,
         0.25854, 0.27001, 0.27667, 0.05146, -0.44994, -0.45381, -0.44362,
         -0.41888, -0.37613, -0.32251, -0.25931, -0.19118, -0.11873, -0.05590, 0]
ref_x = [0, 0.0151, 0.0308, 0.0454, 0.0600, 0.0747, 0.0902, 0.1049, 0.1206,
         0.1352, 0.1450, 0.5005, 0.8501, 0.8647, 0.8804, 0.8950, 0.9106, 0.9253,
         0.9399, 0.9546, 0.9702, 0.9849, 1]

# Problem 2 plot
fig, ax = plt.subplots(1, 2)
fig.set_figwidth(9)
fig.set_figheight(3.5)
for N in Ns:
    ax[0].plot(y[N], u[N], label='{}x{}'.format(N, N), linewidth=1)
    ax[1].plot(x[N], v[N], linewidth=1)

ax[0].plot(ref_y, ref_u, '.k', label='Roy et. al (2015)', markersize=4)
ax[1].plot(ref_x, ref_v, '.k', markersize=4)
ax[0].set_xlabel(r'Normalized $y$')
ax[0].set_ylabel(r'Normalized $u$')
ax[1].set_xlabel(r'Normalized $x$')
ax[1].set_ylabel(r'Normalized $v$')
ax[0].grid()
ax[1].grid()
handles, labels = ax[0].get_legend_handles_labels()
lgd = ax[0].legend(handles, labels, loc='lower center', bbox_to_anchor=(1.0, -0.29),
                      ncol=7, fontsize=9)
fig.tight_layout()
fig.savefig('results/p2.pdf', bbox_inches='tight', bbox_extra_artists=(lgd,))
