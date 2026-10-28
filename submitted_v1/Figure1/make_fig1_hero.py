"""
Fig. 1 (hero, single column, 3 panels):
  (a) phase diagram of the time-averaged boundary QFI in the (h/J, gamma) plane
  (b) QFI dynamics in the topological phase (gamma=0.3, h=0.5J)
  (c) QFI dynamics in the trivial phase   (gamma=0.8, h=1.8J)

Reads:  fig1_*.txt
Writes: Figure1.pdf
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

CMAP = 'RdPu'
DPI = 300
DATA = '.'
OUT = '.'
os.makedirs(OUT, exist_ok=True)

# --- Load data ---
times = np.loadtxt(os.path.join(DATA, 'fig1_times.txt')).flatten()
t_max = times[-1]

h_qfi = np.loadtxt(os.path.join(DATA, 'fig1_phase_diagram_h.txt')).flatten()
g_qfi = np.loadtxt(os.path.join(DATA, 'fig1_phase_diagram_gamma.txt')).flatten()
Z_qfi = np.loadtxt(os.path.join(DATA, 'fig1_phase_diagram_Z.txt'))

data_topo = np.loadtxt(os.path.join(DATA, 'fig1_dynamics_topological.txt'))
data_triv = np.loadtxt(os.path.join(DATA, 'fig1_dynamics_trivial.txt'))
L_dyn = data_topo.shape[1]
vmax_shared = max(data_topo[:, 1:].max(), data_triv[:, 1:].max())

# --- Layout: 3.4 inch single column ---
fig = plt.figure(figsize=(3.4, 4.2))
gs = GridSpec(2, 3, figure=fig,
              height_ratios=[1, 0.85],
              width_ratios=[1, 1, 0.05],
              hspace=0.35, wspace=0.30)

# (a) Phase diagram — top, spanning both columns
ax_a = fig.add_subplot(gs[0, 0:2])
cp = ax_a.pcolormesh(h_qfi, g_qfi, Z_qfi, cmap=CMAP, shading='auto')
cbar_a = fig.colorbar(cp, cax=fig.add_subplot(gs[0, 2]))
cbar_a.set_label(r'$\overline{\mathcal{I}}_1$', fontsize=11)
cbar_a.ax.tick_params(labelsize=7)
ax_a.axvline(1.0, color='black', linestyle='--', lw=1.5, label=r'$|h|=J$')
ax_a.legend(loc='upper right', fontsize=6, facecolor='white',
            edgecolor='none', labelcolor='black', framealpha=0.8)
ax_a.set_xlabel(r'Transverse field ($h/J$)')
ax_a.set_ylabel(r'Anisotropy ($\gamma$)')
ax_a.text(0.14, 0.08, 'Topological', transform=ax_a.transAxes,
          fontsize=6, color='black', ha='center', fontweight='bold')
ax_a.text(0.88, 0.08, 'Trivial', transform=ax_a.transAxes,
          fontsize=6, color='black', ha='center', fontweight='bold')
ax_a.text(0.03, 0.95, '(a)', transform=ax_a.transAxes,
          fontsize=10, color='white', ha='left', va='top', fontweight='bold')

# (b) Topological dynamics — bottom-left
ax_b = fig.add_subplot(gs[1, 0])
im_b = ax_b.imshow(data_topo, aspect='auto', origin='lower',
                   extent=[1, L_dyn, 0, t_max], cmap=CMAP, vmin=0, vmax=vmax_shared)
ax_b.set_xlim(-4, L_dyn)  # margin left of the first site: the locked boundary streak stays visible
ax_b.set_xlabel(r'Lattice site $j$')
ax_b.set_ylabel(r'Time ($tJ$)')
ax_b.text(0.05, 0.92, '(b)', transform=ax_b.transAxes,
          fontsize=10, color='black', ha='left', va='top', fontweight='bold')

# (c) Trivial dynamics — bottom-right
ax_c = fig.add_subplot(gs[1, 1])
im_c = ax_c.imshow(data_triv, aspect='auto', origin='lower',
                   extent=[1, L_dyn, 0, t_max], cmap=CMAP, vmin=0, vmax=vmax_shared)
ax_c.set_xlim(-4, L_dyn)
ax_c.set_xlabel(r'Lattice site $j$')
ax_c.set_yticklabels([])  # share y-axis label with (b)
ax_c.text(0.05, 0.92, '(c)', transform=ax_c.transAxes,
          fontsize=10, color='black', ha='left', va='top', fontweight='bold')

# Shared colorbar for dynamics
cax_bc = fig.add_subplot(gs[1, 2])
cb_bc = fig.colorbar(im_b, cax=cax_bc)
cb_bc.set_label(r'$\mathcal{I}_j$', fontsize=11)
cb_bc.ax.tick_params(labelsize=7)

figname = os.path.join(OUT, 'Figure1.pdf')
plt.savefig(figname, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Saved: {figname}")
