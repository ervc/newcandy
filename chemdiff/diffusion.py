import numpy as np

from . import constants as const
from .utils import grain_abun2dg
from .column import Column

def do_diffusion(col: Column, diffdt: float) -> None:
    diff = col.alpha*col.cs*col.h
    beta = 0.5*diff*diffdt*const.YR
    nzs = col.ncells
    col_abuns = col.get_abundance_array()
    newarray = {}
    for spec in col_abuns:
        newarray[spec] = np.zeros(nzs)
        for j in range(nzs):
            rho_j = col.cells[j].rho
            sp = col_abuns[spec]
            fp = 0
            if j < nzs-1:
                rho_j1 = col.cells[j+1].rho
                fp = 0.5*(rho_j+rho_j1)*beta*(sp[j+1]-sp[j])/col.dz
            elif j == nzs-1:
                fp = 0.
            if j==0:
                fm = -fp
            newarray[spec][j] = sp[j]+(fp-fm)/col.dz/rho_j
            fm = fp
    for spec in newarray:
        for j in range(nzs):
            col.cells[j].abundances[spec] = newarray[spec][j]
            col.cells[j].dust_gas_ratio = grain_abun2dg(newarray['grain'][j],
                                            col.cells[j].grain_size)

def grow_grains(
        col: Column, diffdt: float, pebcomp: dict, 
        growth_timescale_factor: float, growth_height: float) -> dict:
    nzs = col.ncells
    for j in range(nzs):
        cell = col.cells[j]
        t_grow = growth_timescale_factor/(cell.dust_gas_ratio*col.omega)
        if cell.z/col.h <= growth_height:
            deps = -diffdt*const.YR/t_grow
        else:
            deps = 0
        # if spec is an ice then reduce its abundance by deps
        # and add that missing ice to the pebcomp dictionary
        for spec in cell.abundances:
            if spec.startswith('grain'):
                d_ice = deps*cell.abundances[spec]
                cell.abundances[spec] += d_ice
                if spec not in pebcomp:
                    pebcomp[spec] = 0
                pebcomp[spec] -= d_ice*cell.nh*col.dz
        cell.dust_gas_ratio = grain_abun2dg(cell.abundances['grain'],
                                            col.cells[j].grain_size)
    return pebcomp

    