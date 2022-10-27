import numpy as np

from . import constants as const

def get_omega(r: float) -> float:
    return np.sqrt(const.G*const.MSUN/r)/r

def get_midplane_temp(r: float) -> float:
    # from Krijt et al. 2018
    return 130*(r/const.AU)**(-1/2)

def get_soundspeed(r: float) -> float:
    tmid = get_midplane_temp(r)
    return np.sqrt(const.BK*tmid/const.MBAR)

def get_scaleheight(r: float) -> float:
    om = get_omega(r)
    cs = get_soundspeed(r)
    return cs/om

def get_surface_density(r: float) -> float:
    Mdisk=0.05*const.MSUN
    p = 1
    rc = 100*const.AU
    sigc = (2-p)*Mdisk/(2*np.pi*rc*rc)
    sig = sigc*(r/rc)**(-p) * np.exp(-(r/rc)**(2-p))
    return sig

def get_temperature(r: float, z: float) -> float:
    tmid = get_midplane_temp(r)
    tatm = 2*tmid
    h = get_scaleheight(r)
    zq = 3 # scale height at which temp=tatm
    temp = tmid+(tatm-tmid)*(np.sin(np.pi*z/2./zq/h))**(4)
    if z > zq*h:
        temp = tatm
    return temp

def get_density(r: float, z: float) -> float:
    sig = get_surface_density(r)
    h = get_scaleheight(r)
    rho0 = sig/np.sqrt(2.*np.pi)/h
    rho = rho0*np.exp(-z*z/2./h/h)
    return rho
