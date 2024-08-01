import sys
from .constants import *

class Cell(object):
    '''
    A cell has:

    a physical location
        r,z

    physical parameters for the solver
        chi, cosmic, grain size, dust/gas ratio

    source model parameters
        given
            Av, nh, Tgas, Tdust, xrays
        calculated
            NCO, NH2, (density)

    abundances
        dict of abundances

    diffusion parameters
        given 
            omega, 
            alpha, dz, dt (solver params?)
        calculated
            cs, h, D, beta
    '''


    def __init__(self,r,z,
            chi=1,cosmic=1.3e-17,grain_size=0.1,dust_gas_ratio=0.01,
            av=1,rho=1e10,Tgas=50,Tdust=50,xray=0,NCO=1.,NH2=1.,NHD=1.,NH=1.,
            abundances = {}):
        self.r = r # cm
        self.z = z # cm

        # phys
        self.chi = chi
        self.cosmic = cosmic # s-1
        self.grain_size = grain_size # micron
        self.dust_gas_ratio = dust_gas_ratio

        self.av = av # mag
        self.rho = rho # g.cm-3
        self.nh = 2*rho/MBAR # cm-3
        self.Tgas = Tgas # K
        self.Tdust = Tdust # K
        self.xray = xray # s-1
        self.NCO = NCO # cm-2
        self.NH2 = NH2 # cm-2
        self.NHD = NHD # cm-2
        self.NH = NH # cm-2

        self.abundances = dict(abundances)


    def write_chem_inputs(self,tf: float,abs_err: float,rel_err: float,abun_out: str='all',
            f_net: str='network.chm',f_input: str='input.ini',f_source: str='source.mdl'):
        """Write the chmical input and source files for astrochem

        Args:
            tf (float): final time for chemistry
            abs_err (float): absolute error for chemistry. If abundance is less than this set to zero
            rel_err (float): relative error for chemistry
            abun_out (str, optional): what species to track. Defaults to 'all'.
            f_net (str, optional): chemical network to use for chemistry. Defaults to 'network.chm'.
            f_input (str, optional): name of input file to write. Defaults to 'input.ini'.
            f_source (str, optional): name of source file to write. Defaults to 'source.mdl'.
        """
        
        with open(f_input,'w') as f:
            f.write('[files]\n')
            f.write(f'source = {f_source}\n'+
                f'chem = {f_net}\n')

            f.write('[phys]\n')
            f.write(f'chi = {self.chi:.3e}\n'+
                f'cosmic = {self.cosmic:.2e}\n'+
                f'grain_size = {self.grain_size:.2e}\n'+
                f'grain_gas_mass_ratio = {self.dust_gas_ratio:.2e}\n')

            f.write('[solver]\n')
            f.write('ti = 1.00e-06\n'+
                f'tf = {tf:.2e}\n'+
                f'abs_err = {abs_err:.1e}\n'+
                f'rel_err = {rel_err:.1e}\n')

            f.write('[abundances]\n')
            for spec in self.abundances:
                try:
                    if self.abundances[spec] > abs_err:
                        f.write(f'{spec} = {self.abundances[spec]:.15e}\n')
                except:
                    print('Problem with Abundances:\n',self.abundances)
                    sys.exit()


            f.write('[output]\n')
            f.write(f'abundances = {abun_out}\n'+
                'time_steps = 64\n'+
                'trace_routes = 0')

        # print(f'Input written to {f_input}')

        with open(f_source,'w') as f:
            f.write('# self Av[mag] n(H)[cm-3] Tgas[K] Tdust[K] NCO[cm-2] NH2[cm-2] NHD[cm-2] xray-ion[s-1] R[au] Z[au]\n')
            f.write(f'0    {self.av:.3e}    {self.nh:.3e}    {self.Tgas:.3e}    {self.Tdust:.3e}    {self.NCO:.3e}    {self.NH2:.3e}    {self.NHD:.3e}    {self.xray:.3e}    {self.r/AU:.2f}    {self.z/AU:.2f}')

        # print(f'Source written to {f_source}')

    def update_abundances(self,new_dict: dict):
        """Update the abundances of each species in the cell

        Args:
            new_dict (dict): dictionary of new species and abundances

        """
        self.abundances=dict(new_dict)