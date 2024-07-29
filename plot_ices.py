import numpy as np
import matplotlib.pyplot as plt

import chemdiff as cd

def get_pebble_abuns(cdout, pebfile):
    times = cdout['times']
    numdens = cdout['density'] # number density of hydrogen cm-3
    heights = cdout['heights']
    nz = len(heights)
    dz = heights[1]-heights[0] # dz is constant

    totNH = 0
    for z in range(nz):
        totNH += numdens[z]*dz

    pebcomps = np.zeros_like(times)
    with open(pebfile) as f:
        t=0
        for i,line in enumerate(f):
            if i==0: continue # skip first line
            time,ice,abun = line.split()
            if ice == 'grainH2O':
                pebcomps[t] = float(abun)*totNH
                t+=1

    return pebcomps

def get_ice_abuns(cdout):
    # candy ouputs
    specs = cdout['species']
    abuns = cdout['abundances']
    numdens = cdout['density'] # number density of hydrogen cm-3
    heights = cdout['heights']
    nz = len(heights)
    dz = heights[1]-heights[0] # dz is constant

    spec = 'grainH2O'
    subtotal = 0
    k = np.where(specs==spec)[0][0]
    for z in range(nz):
        subtotal += abuns[:,z,k]*numdens[z]*dz

    return subtotal

def get_alt_ice_abuns(cdout):
    specs = cdout['species']
    abuns = cdout['abundances']
    numdens = cdout['density'] # number density of hydrogen cm-3
    heights = cdout['heights']
    nz = len(heights)
    dz = heights[1]-heights[0] # dz is constant

    spec = 'grainH2O'
    subtotal = 0
    k = np.where(specs==spec)[0][0]
    kdust = np.where(specs=='grain')[0][0]
    for z in range(nz):
        subtotal += abuns[:,z,k]*numdens[z]*dz - 

    return subtotal


def get_gas_abuns(cdout):
    # candy ouputs
    specs = cdout['species']
    abuns = cdout['abundances']
    numdens = cdout['density'] # number density of hydrogen cm-3
    heights = cdout['heights']
    nz = len(heights)
    dz = heights[1]-heights[0] # dz is constant

    spec = 'H2O'
    subtotal = 0
    k = np.where(specs==spec)[0][0]
    for z in range(nz):
        subtotal += abuns[:,z,k]*numdens[z]*dz

    return subtotal


def get_analytic_start(cdout):
    iceabun = 2.e-5
    numdens = cdout['density'] # number density of hydrogen cm-3
    heights = cdout['heights']
    nz = len(heights)
    dz = heights[1]-heights[0] # dz is constant

    tot = 0.
    for z in range(nz):
        tot += iceabun*numdens[z]*dz

    return tot

def read_astrochem(hfile):
    return cd.candyio.get_abundict(hfile)

def plot_one_static(hfile,outfile):
    times,abundict = read_astrochem(hfile)
    fig,axs = plt.subplots(2,1,sharex=True)
    ax=axs[0]
    ax.plot(times,abundict['H2O'],ls='-')
    ax.plot(times,abundict['grainH2O'],ls='--')
    total = abundict['H2O']+abundict['grainH2O']
    ax.plot(times,total,ls=':')
    ax.set(xlabel='time',ylabel='abundance',xscale='log',yscale='log')

    ax=axs[1]
    ax.plot(times,total/total[0])
    ax.set(xlabel='time',ylabel='total/total_0')

    plt.savefig(outfile,bbox_inches='tight')
    return

def plot_all_static(cddir,outfile):
    candyfile = cddir+'/output_abuns.npz'
    cdout = np.load(candyfile)
    heights = cdout['heights']
    nz = len(heights)
    dz = heights[1]-heights[0] # dz is constant

    totalice = 0
    totalgas = 0
    for z in range(nz):
        hfile   = cddir+f'/z{z:0>2}/astrochem_output.h5'
        srcfile = cddir+f'/z{z:0>2}/source.mdl'
        nh,_ = cd.candyio.get_nh_and_av(srcfile)
        times,abuns = read_astrochem(hfile)

        
        totalgas += abuns['H2O']*nh*dz
        totalice += abuns['grainH2O']*nh*dz
    totaltot = totalgas+totalice

    fig,axs = plt.subplots(2,1,sharex=True)
    ax=axs[0]
    ax.plot(times,totalgas,ls='-')
    ax.plot(times,totalice,ls='--')
    ax.plot(times,totaltot,ls=':')
    ax.set(ylabel='col density',xscale='log',yscale='log')

    ax=axs[1]
    ax.plot(times,totaltot/totaltot[0])
    ax.set(xlabel='time',ylabel='total/total_0')

    plt.savefig(outfile,bbox_inches='tight')
    return

def main():
    # candy output file
    outdir = 'examples/ice_check/'
    
    # plot_one_static(outdir+'/z00/astrochem_output.h5',outdir+'z0_abuns.png')
    # plot_all_static(outdir,outdir+'/allstatic.png')
    # return

    # load candy results
    

    # candy ouputs
    
    fig,axs = plt.subplots(2,1,sharex=True,height_ratios=(4,1))

    outdirs = [
        'examples/ice_check/',
        'examples/ice_check_difffirst',
    ]

    lss = ['-','--']

    for ls,outdir in zip(lss,outdirs):

        candyfile = outdir+'/output_abuns.npz'
        pebfile = outdir+'/pebcomp.out'

        cdout = np.load(candyfile)

        times = cdout['times']


        gasabun = get_gas_abuns(cdout)
        iceabun = get_ice_abuns(cdout)
        pebabun = get_pebble_abuns(cdout,pebfile)

        totabun = gasabun+iceabun+pebabun

        analytic_abun = get_analytic_start(cdout)

        ax=axs[0]
        # ax.plot(times,gasabun,label='gas',ls=ls)
        ax.plot(times,iceabun,label='ice',ls=ls)
        ax.plot(times,pebabun,label='peb',ls=ls)
        ax.plot(times,totabun,label='total',ls=ls)

        print('max value, min value = ',np.max(totabun),np.min(totabun))

        ax.legend(loc='upper left',bbox_to_anchor=(1,1))
        ax.set(xlabel='Time',ylabel='Column total number surface density',
            xscale='log',yscale='log')
        
        ax=axs[1]
        ax.axhline(1,ls=':',c='k')
        ax.plot(times,totabun/analytic_abun,label=outdir,ls=ls)
        ax.set(xlabel='Time',ylabel='Total/Total_0')
        ax.legend()


    plt.savefig('examples/ice_check/ices_compare.png',bbox_inches='tight')

if __name__ == '__main__':
    main()