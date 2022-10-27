import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import chemdiff as cd

candyfile = 'static_example/output_abuns.npz'
anifile = 'static_example/animation.gif'

cdout = np.load(candyfile)

plotspecs = ['CO','grainCO','H2O','grainH2O','grain']
specs = cdout['species']
abuns = cdout['abundances']
times = cdout['times']

r = 30*cd.constants.AU
h = cd.disk.get_scaleheight(r)

fig,ax = plt.subplots()
lines = {}
for spec in plotspecs:
    if spec in specs:
        k = np.where(specs==spec)[0][0]
        line, = ax.plot(abuns[0,:,k],cdout['heights']/h,label=spec)
        lines[spec] = line
ax.set(xscale='log',xlabel='abun',ylabel='Z/H')
ax.legend()
title = ax.set_title(f'{times[0]:.0f} yr')

def animate(i):
    for spec,line in lines.items():
        k = np.where(specs==spec)[0][0]
        line.set_xdata(abuns[i,:,k])
    title.set_text(f'{times[i]:.0f} yr')

anim = FuncAnimation(fig,animate,frames=len(times))
anim.save(anifile)