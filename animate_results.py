# standard imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# local imports
import chemdiff as cd

# candy output file
candyfile = 'growth_example/output_abuns.npz'
# the name of the gif we will be making
anifile = 'growth_example/animation.gif'

# load candy results
cdout = np.load(candyfile)

# which species do we want to plot
plotspecs = ['CO','grainCO','H2O','grainH2O','grain']

# candy ouputs
specs = cdout['species']
abuns = cdout['abundances']
times = cdout['times']

# some helpful physical parameters
r = 30*cd.constants.AU
h = cd.disk.get_scaleheight(r)

# create the figure
fig,ax = plt.subplots()
# create an empty dictionary to store the lines that are plotted
# the dictionary will have entries so that the line for H2O
# can be indexed using line[H2O]
lines = {}
# loop through, get the abundances for the species we want to plot
# then plot the initial abundances and store the line created
# in the dictionary
for spec in plotspecs:
    if spec in specs:
        # find which species index corresponds with the species of interest
        k = np.where(specs==spec)[0][0]
        # plot abundance at time 0, all heights, species k vs Z/H
        # ax.plot returns a tuple, but we only care about the line object,
        # so that's why I say line, = ...
        line, = ax.plot(abuns[0,:,k],cdout['heights']/h,label=spec)
        lines[spec] = line
# axes labels, title and legend
ax.set(xscale='log',xlabel='abun',ylabel='Z/H')
ax.legend()
title = ax.set_title(f'{times[0]:.0f} yr')

# the animate function is called when the animation is made, where the
# argument i is the frame number
def animate(i):
    # loop through the species and lines from the lines dictionary
    for spec,line in lines.items():
        # Again, get the right index
        k = np.where(specs==spec)[0][0]
        # adjust the x data (since I'm plotting height as y, which doesn't change)
        # the new x data is abundance at time i, all heights, species index k
        line.set_xdata(abuns[i,:,k])
    # adjust the title
    title.set_text(f'{times[i]:.0f} yr')
    # Technically you're supposed to return the objects you modify, I think
    # that helps it run faster, but it works fine if you just don't return
    # anything here

# create and save the animation, I'm passing frames=len(times) to assure
# we loop through all the times from the output
anim = FuncAnimation(fig,animate,frames=len(times))
anim.save(anifile)