"""
Fig. 6 (single column): interaction robustness.
Time-averaged boundary QFI vs interaction strength delta/J from exact
statevector simulations, for (gamma, h) = (1, 0), (1, 0.2J), (0.5, 0),
(0.5, 0.2J) and system sizes L = 12, 16, 20.
R_y encoding, k=j=1, averaging window tJ in [25, 50].

Reads:  fig6_interacting_L*_g*_h*.txt
Writes: Figure6.pdf
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.size': 8,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

DPI = 300
DATA = '.'
OUT = '.'
os.makedirs(OUT, exist_ok=True)

param_sets = [(1.0, 0.0), (1.0, 0.2), (0.5, 0.0), (0.5, 0.2)]
L_list = [12, 16, 20]

# same family as Figure 5: pinks for gamma=1, purples for gamma=0.5,
# dark for h=0, light for h=0.2J
palette = {
    (1.0, 0.0): '#C71585',
    (1.0, 0.2): '#FF80C0',
    (0.5, 0.0): '#4B0082',
    (0.5, 0.2): '#9B6FD0',
}
param_labels = {
    (1.0, 0.0): r'$\gamma{=}1,\; h{=}0$',
    (1.0, 0.2): r'$\gamma{=}1,\; h{=}0.2J$',
    (0.5, 0.0): r'$\gamma{=}0.5,\; h{=}0$',
    (0.5, 0.2): r'$\gamma{=}0.5,\; h{=}0.2J$',
}
L_styles = {12: (':', 'o', 0.45), 16: ('--', 's', 0.65), 20: ('-', 'o', 1.0)}

fig, ax = plt.subplots(figsize=(3.4, 2.6))
fig.subplots_adjust(left=0.14, right=0.97, bottom=0.17, top=0.95)

for (gamma, h_field) in param_sets:
    color = palette[(gamma, h_field)]
    lbl = param_labels[(gamma, h_field)]
    for N in L_list:
        fname = os.path.join(
            DATA, f'fig6_interacting_L{N}_g{gamma:.2f}_h{h_field:.2f}.txt')
        if not os.path.exists(fname):
            continue
        dat = np.loadtxt(fname)
        delta, fq = dat[:, 0], dat[:, 1]
        ls, mk, alpha = L_styles[N]
        label = lbl if N == L_list[-1] else None
        ax.plot(delta, fq, ls=ls, marker=mk, color=color,
                alpha=alpha, lw=1.5, markersize=3,
                markeredgewidth=0.4, markeredgecolor='k',
                label=label)

# L legend entries (gray)
for N in L_list:
    ls, mk, alpha = L_styles[N]
    ax.plot([], [], ls=ls, marker=mk, color='gray', alpha=alpha,
            lw=1.5, markersize=3, markeredgewidth=0.4, markeredgecolor='k',
            label=f'$L={N}$')

ax.set_xlabel(r'Interaction Strength ($\delta / J$)')
ax.set_ylabel(r'$\overline{\mathcal{I}}_1$')
ax.set_xlim(0, 0.8)
ax.set_ylim(0, 1.05)
ax.legend(loc='lower left', fontsize=6.5, ncol=2,
          columnspacing=0.6, handlelength=1.5,
          framealpha=0.9, edgecolor='none')

figname = os.path.join(OUT, 'Figure6.pdf')
plt.savefig(figname, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Saved: {figname}")
