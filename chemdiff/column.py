import numpy as np

from . import disk

class Column:
    """1D Column of cells to allow for chemistry, diffusion, and grain growth
    """
    def __init__(self,r: float,alpha: float=1e-3,ncells: int=50):
        """Column object for chemistry and diffusion within CANDY.
        Contains 1d (vertical) array of cells at a given location in the disk.

        Args:
            r (float): radial location of the column [cm]
            alpha (float, optional): alpha viscosity param. Defaults to 1e-3.
            ncells (int, optional): number of cells to divide column. Defaults to 50.
        """
        self.r = r # cm
        self.ncells = ncells
        self.tmid = disk.get_midplane_temp(r) # K
        self.omega = disk.get_omega(r) #1/s
        self.cs = disk.get_soundspeed(r) #cm/s
        self.h = disk.get_scaleheight(r) # cm

        self.alpha = alpha
        self.dz = 5*self.h/ncells # cm

        self.cells = np.empty(ncells,dtype='object')

    def get_abundance_array(self) -> dict:
        """Returns dictionary of chemical abundances for each cell

        Returns:
            dict: dictionary of abundances. Keys are the string of the chemical and values are arrays of length (nz,)
        """
        all_abunds = {}
        for i in range(self.ncells):
            cell = self.cells[i]
            for spec in cell.abundances:
                if spec not in all_abunds:
                    all_abunds[spec] = np.zeros((self.ncells))
                all_abunds[spec][i] = cell.abundances[spec]
        self.all_abunds = all_abunds
        return all_abunds

    def update_column_densities(self,opacity: float):
        """Update the column densities and visual extinctions throughout
        the column given the current cell abundances. This should be
        called after column abundances are changed to remain consistent

        Args:
            opacity (float): The opacity of small dust grains [cm2 g-1]ß
        """
        NCO = 0.
        NH2 = 0.
        NHD = 0.
        NH = 0.
        tau = 0.
        xray = 0.
        for j in reversed(range(self.ncells)):
            cell = self.cells[j]
            rho = cell.rho
            nh = cell.nh
            nco = 0.
            if 'CO' in cell.abundances:
                nco = nh*cell.abundances['CO']
            nh2 = 0.
            if 'H2' in cell.abundances:
                nh2 = nh*cell.abundances['H2']
            nhd = 0
            if 'HD' in cell.abundances:
                nhd = cell.abundances['HD']*nh
            NCO += nco*self.dz
            NH2 += nh2*self.dz
            NHD += nhd*self.dz
            NH += nh*self.dz
            tau += rho*opacity*self.dz*cell.dust_gas_ratio
            cell.NCO = NCO
            cell.NH2 = NH2
            cell.NHD = NHD
            cell.NH = NH
            cell.av = tau/3.02