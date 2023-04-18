import numpy as np
from . import constants

def grain_abun2dg(grain_abun: float,grain_size=0.1,grain_density=1.81) -> float:
    """Convert grain abundance to a dust-to-gas mass ratio"""
    grain_size = grain_size*1e-4 # convert micron to cm
    dg = grain_abun*grain_density*(4./3.)*np.pi*pow(grain_size,3.)/constants.MH
    return dg

def dg2grain_abun(dg: float,grain_size=0.1,grain_density=1.81) -> float:
    """Convert dust-to-gas mass ratio to a grain abundance"""
    grain_size = grain_size*1e-4 # convert micron to cm
    grain_abun = dg*constants.MH/(grain_density*(4./3.)*np.pi*pow(grain_size,3.))
    return grain_abun