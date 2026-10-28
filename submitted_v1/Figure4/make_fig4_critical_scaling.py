"""
Fig. 4 (single column): critical scaling of the boundary QFI near the
topological transition. Log-log plot of the time-averaged boundary QFI vs
the distance to criticality (1 - h/J) for L = 50, 100, 200, 400, together
with the asymptotic (1-h/J)^2 guide line. gamma=1.0, k=j=1.

The x axis carries a dense set of labeled log10 ticks
(0.02, 0.05, 0.1, 0.2, 0.5) plus unlabeled logarithmic minor ticks.

Reads:  fig4_critical_scaling.txt
Writes: Figure4.pdf
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, LogLocator, FuncFormatter, NullFormatter
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

d = np.loadtxt(os.path.join(DATA, 'fig4_critical_scaling.txt'))
dist = 1.0 - d[:, 0]           # distance to criticality 1 - h/J
L_list = [50, 100, 200, 400]
colors = ['#FF66B8', '#E6007A', '#A30059', '#4A0033']

fig, ax = plt.subplots(figsize=(3.4, 2.6))
fig.subplots_adjust(left=0.15, right=0.97, bottom=0.17, top=0.95)

for i, (L, c) in enumerate(zip(L_list, colors)):
    ax.loglog(dist, d[:, i + 1], 'o-', color=c, lw=1.2,
              markersize=3, label=f'$L={L}$')

dg = np.linspace(0.02, 0.25, 10)
ax.loglog(dg, 3.8 * dg**2, '--', color='red', lw=1.8,
          label=r'$\propto (1-h/J)^2$')

# Labeled log10 ticks on the x axis + unlabeled logarithmic minor ticks
xticks = [0.02, 0.05, 0.1, 0.2, 0.5]
ax.xaxis.set_major_locator(FixedLocator(xticks))
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1,
                                      numticks=30))
ax.xaxis.set_minor_formatter(NullFormatter())

ax.set_xlabel(r'Distance to Criticality ($1 - h/J$)')
ax.set_ylabel(r'Time-Averaged Boundary QFI ($\overline{\mathcal{I}}_1$)')
ax.legend(loc='lower right', fontsize=7, framealpha=0.9, edgecolor='none')

figname = os.path.join(OUT, 'Figure4.pdf')
plt.savefig(figname, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Saved: {figname}")
