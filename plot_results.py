# import (mostly) standard libraries
import numpy as np
import matplotlib.pyplot as plt

# import chemdiff locally
import chemdiff as cd

# get the names of the CANDY output file and the image we want to make
cdfile = 'examples/growth_example/output_abuns.npz'
pngfile = 'examples/growth_example/results.png'

# load the output file into NpzFile class
cdout = np.load(cdfile)

# we can get the outputs that we care about by indexing 
# like a python dictionary.
# You can see what is indices are available by printing
### in : print(cdout.files)
### out: ['times', 'species', 'abundances', 'heights', 'density', 'extinctions']
specs = cdout['species']    # list of all species in the output file
abuns = cdout['abundances'] # array of all the abundances of each species
times = cdout['times']      # list of all the output times
# abundances has shape (nt,nz,nspec)
# species has shape (nspec)
# times has shape (nt)
# heights has shape (nz)
# density has shape (nz)
# extinctions has shape (nt,nz)

# A list of species we want to plot
plotspecs = ['CO','grainCO','H2O','grainH2O','grain']

# some useful constants to have for plotting
r = 30*cd.constants.AU
# here we're using one of the helpful functions from chemdiff.disk
# these disk parameters can be changed by editing funcitons in chemdiff/disk.py
h = cd.disk.get_scaleheight(r)

# create a figure and axis dictionary using handy subplot_mosaic function
fig,axd = plt.subplot_mosaic(
    """
    ppaa
    ppww
    ppmm
    """, tight_layout=True
)
# Here, I want to plot the profile of a few species of interest at the end of
# simulation.

# for each species we want to plot
for spec in plotspecs:
    # if that species is included in the chemdiff output
    if spec in specs:
        # get the index of the species
        k = np.where(specs==spec)[0][0]
        # plot the output for the last time, all heights, species of interest
        # abundances are indexed as abuns[nt,nz,nspec]
        axd['p'].plot(abuns[-1,:,k],cdout['heights']/h,label=spec)
# format the axes
axd['p'].set(xscale='log',xlabel='abun',ylabel='Z/H')
axd['p'].legend()
title = axd['p'].set_title(f'{times[-1]:.0f} yr')

# Here I want to plot the time evolution of the midplane, warm molecular region,
# and atmosphere of the disk. 
# First I define the three regions of the disk (arbitrarily) as the lowest cell,
# a cell half way up the column, and the cell at the very top
cells = {
    'm' : 0,
    'w' : len(cdout['heights'])//2,
    'a' : len(cdout['heights'])-1,
}

# loop through each of the three cells defined above
for label,cell in cells.items():
    # select which axis we want to draw on
    ax = axd[label]
    # for each species
    for spec in plotspecs:
        # if that species exists
        if spec in specs:
            # get the index
            k = np.where(specs==spec)[0][0]
            # plot all times, the appropriate height, and the right species
            ax.plot(times,abuns[:,cell,k],label=spec)
    # format each axis
    ax.set(xscale='log',ylabel='abun',yscale='log')
    ax.set(title=f'z/h = {cdout["heights"][cell]/h:.2f}')
    if label=='m':
        ax.set(xlabel='time [yr]')
    else:
        ax.set(xticks=[],xticklabels=[])

# save the output!
plt.savefig(pngfile)