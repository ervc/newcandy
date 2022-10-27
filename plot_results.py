import numpy as np
import matplotlib.pyplot as plt

import chemdiff as cd

cdfile = 'diff_example/output_abuns.npz'
pngfile = 'diff_example/results.png'


cdout = np.load(cdfile)

plotspecs = ['CO','grainCO','H2O','grainH2O','grain']
specs = cdout['species']
abuns = cdout['abundances']
times = cdout['times']

r = 30*cd.constants.AU
h = cd.disk.get_scaleheight(r)

fig,axd = plt.subplot_mosaic(
    """
    ppaa
    ppww
    ppmm
    """, tight_layout=True
)
for spec in plotspecs:
    if spec in specs:
        k = np.where(specs==spec)[0][0]
        axd['p'].plot(abuns[-1,:,k],cdout['heights']/h,label=spec)
axd['p'].set(xscale='log',xlabel='abun',ylabel='Z/H')
axd['p'].legend()
title = axd['p'].set_title(f'{times[-1]:.0f} yr')

cells = {
    'm' : 0,
    'w' : len(cdout['heights'])//2,
    'a' : len(cdout['heights'])-1,
}

for label,cell in cells.items():
    ax = axd[label]
    for spec in plotspecs:
        if spec in specs:
            k = np.where(specs==spec)[0][0]
            ax.plot(times,abuns[:,cell,k],label=spec)
    ax.set(xscale='log',ylabel='abun',yscale='log')
    ax.set(title=f'z/h = {cdout["heights"][cell]/h:.2f}')
    if label=='m':
        ax.set(xlabel='time [yr]')
    else:
        ax.set(xticks=[],xticklabels=[])

plt.savefig(pngfile)