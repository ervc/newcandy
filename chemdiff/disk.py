import numpy as np

from . import constants as const

"""Analytic equations for disk structure
"""

def get_omega(r: float) -> float:
    """Return Keplerian frequency at given radial location

    Args:
        r (float): radial distance [cm]

    Returns:
        float: Keplerian Frequency [s-1]
    """
    return np.sqrt(const.G*const.MSUN/r)/r

def get_midplane_temp(r: float) -> float:
    """Midplane temperature of the disk, following Krijt et al. 2018

    Args:
        r (float): radial distance [cm]

    Returns:
        float: midplane temperature [K]
    """
    return 130*(r/const.AU)**(-1/2)

def get_soundspeed(r: float) -> float:
    """Return the soundspeed at a given location in the disk

    Args:
        r (float): radial distance [cm]

    Returns:
        float: soundspeed [cm/s]
    """
    tmid = get_midplane_temp(r)
    return np.sqrt(const.BK*tmid/const.MBAR)

def get_scaleheight(r: float) -> float:
    """Return the scaleheight at a given location

    Args:
        r (float): radial distance [cm]

    Returns:
        float: scaleheight [cm]
    """
    om = get_omega(r)
    cs = get_soundspeed(r)
    return cs/om

def get_surface_density(r: float) -> float:
    """Surface density at a given location in the disk following a powerlaw. Assumes disk is 0.05Msun in mass, and edge of disk is at 100 AU. Power law slope of -1

    Args:
        r (float): radial distance [cm]

    Returns:
        float: Surface density, Sigma(r) [g cm-2]
    """
    Mdisk=0.05*const.MSUN
    p = 1
    rc = 100*const.AU
    sigc = (2-p)*Mdisk/(2*np.pi*rc*rc)
    sig = sigc*(r/rc)**(-p) * np.exp(-(r/rc)**(2-p))
    return sig

def get_temperature(r: float, z: float) -> float:
    """Get the temperature at a given 2d location for an irradiated disk

    Args:
        r (float): radial distance [cm]
        z (float): height above midplane [cm]

    Returns:
        float: Disk Temperature [K]
    """
    tmid = get_midplane_temp(r)
    tatm = 2*tmid
    h = get_scaleheight(r)
    zq = 3 # scale height at which temp=tatm
    temp = tmid+(tatm-tmid)*(np.sin(np.pi*z/2./zq/h))**(4)
    if z > zq*h:
        temp = tatm
    return temp

def get_density(r: float, z: float) -> float:
    """Return the gas density at a given 2d location

    Args:
        r (float): radial distance [cm]
        z (float): height above midplane [cm]

    Returns:
        float: gas density [g cm-3]
    """
    sig = get_surface_density(r)
    h = get_scaleheight(r)
    rho0 = sig/np.sqrt(2.*np.pi)/h
    rho = rho0*np.exp(-z*z/2./h/h)
    return rho
