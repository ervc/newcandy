import configparser
import subprocess

import numpy as np
import h5py as h5

from .column import Column

def get_touts(tf: float) -> list[float]:
    """Create the default list of output times

    Args:
        tf (float): final time to create touts until

    Returns:
        list[float]: list of output times
    """
    touts = []
    base = 1.0
    exp = 3
    time = 0
    while time < tf:
        time = base*10**exp
        touts.append(time)
        
        if exp < 5:
            base += 1.0
        else:
            base += 0.5

        if base >= 10:
            base = base/10
            exp += 1
    return touts

def parse_list(s: str) -> list:
    """Read in a string and convert to a python list. Used to read in
    touts from input file.

    Args:
        s (str): list string

    Returns:
        list: parsed list
    """
    assert s.startswith('[') and s.endswith(']'), "touts format must be as a list with square brackets"
    ts = list(map(float,s[1:-1].split(',')))
    return ts


def get_defaults() -> tuple[dict, dict, dict]:
    """Returns the default parameters for CANDY run as a tuple of
    dictionaries.

    Returns:
        tuple[dict, dict, dict]: tuple of dictionaries containing default model params, physical params, and abundances
    """
    # default model parameters
    model_inputs = {
        'chmfile' : 'astrochem/networks/umist12_x.chm',
        'pebfile' : 'pebcomp.out',
        'ti' : 0,
        'tf' : 1e6,
        'chemtime' : 500.0,
        'difftime' : 100.0,
        'outputdir' : 'r00',
        'outfile' : 'output_abuns.npz',
        'ncells' : 50
    }
    model_inputs['touts'] = get_touts(model_inputs['tf'])
    # default physical paramters
    phys_inputs = {
        'r' : 30,
        'r_units' : 'au',
        'alpha' : 1.e-3,
        'chi' : 50.,
        'cosmic' : 1.3e-17,
        'grain_size' : 0.1,
        'dg0' : 0.01,
        'opacity' : 1.0e5,
        'growth_timescale_factor' : 1.,
        'growth_height' : 1.,
        'growth_delay_time' : 0.,
    }
    # default initial chemical abundances
    abundances = {
        'H2' : 0.5,
        'He' : 9.75e-2,
        'NH3': 1.45e-6,
        'H2O': 1.18e-4,
        'CO' : 6.00e-5,
        'N2' : 2.00e-5,
        'CH4': 2.00e-6,
        'CH3OH' : 1.00e-6,
        'H2S' : 1.91e-8,
        'CO2' : 5.00e-5,
        'HCN' : 3.50e-7,
        'grain' : 2.2e-12
    }
    return model_inputs,phys_inputs,abundances

def read_infile(fin: str) -> tuple[dict, dict, dict]:
    """Reads the input file to create input parameters dictionaries for CANDY

    Args:
        fin (str): name of input file

    Returns:
        tuple[dict, dict, dict]: dictionaries of model params, physical params, and abundances
    """
    config = configparser.ConfigParser()
    # the next line retains case sesnitivity for input parameters
    config.optionxform = str
    config.read(fin)
    model,phys,abuns = get_defaults()
    for key in config['model']:
        if key in ['chmfile','pebfile','outputdir','outfile']:
            model[key] = config['model'][key]
        elif key in ['touts']:
            model[key] = parse_list(config['model'][key])
        elif key in ['ncells']:
            model[key] = int(config['model'][key])
        else:
            model[key] = float(config['model'][key])
    if 'touts' not in config['model'].keys():
        model['touts'] = get_touts(model['tf'])
    if model['tf'] not in model['touts']:
        print('Warning: tf is not in touts, including tf in list of output times')
        model['touts'].append(model['tf'])

    dterror = 'Diffusion time must be less than or equal to chemtime.\n' \
            +f'chemtime = {model["chemtime"]}\n' \
            +f'difftime = {model["difftime"]}'
    assert model['chemtime'] >= model['difftime'], dterror
    
    model['touts'].sort()
    print('touts = ',model['touts'])
    for key in config['phys']:
        if key != 'r_units':
            phys[key] = float(config['phys'][key])
        else:
            phys[key] = config['phys'][key]
    if 'abundances' in config:
        abuns = {key:float(val) for (key,val) in config['abundances'].items()}
    return model,phys,abuns

def save_outputs(col: Column, nout: int, outputdirr: str) -> None:
    """Copy and paste astrochem outputs to store as backup or to continue a run

    Args:
        col (Column): CANDY column of cells
        nout (int): the output index
        outputdirr (str): directory where savefiles should be stored
    """
    for j in range(col.ncells):
        dirr = f'{outputdirr}/z{j:0>2}'
        subprocess.run(['cp','astrochem_output.h5',
            f'astrochem_output_t{nout}.h5'], cwd=dirr)
        subprocess.run(['cp','source.mdl',f'source_t{nout}.mdl'],cwd=dirr)
        subprocess.run(['cp','input.ini',f'input_t{nout}.ini'],cwd=dirr)

def new_save_outputs(col: Column, nout: int, time: float,
                     outputdirr: str, outputfile: str,
                     delete_subdir: bool = False) -> None:
    """Save the current abundances of ice and gas into npz file. This
    saves from the memory and accounts for ices lost due to grain growth

    Args:
        col (Column): CANDY column of cells
        nout (int): output index
        time (float): the current time
        outputdirr (str): output directory
        outputfile (str): output file name
        delete_subdir (bool, optional): Delete subdirectories after values are saved. Defaults to False.
    """
    
    save_outputs(col,nout,outputdirr)
    nzs = col.ncells
    if nout==0:
        species = get_specs_array(outputdirr)
        nspec = len(species)
        times = read_touts_from_file(f'{outputdirr}/times.out')
        nts = len(times)

        gas_dens = np.zeros(nzs)
        extinctions = np.zeros((nts,nzs))
        abundances = np.zeros((nts,nzs,nspec))
        heights = np.zeros(nzs)

        for j in range(nzs):
            cell = col.cells[j]
            gas_dens[j] = cell.nh
            extinctions[0,j] = cell.av
            heights[j] = cell.z
            for k in range(nspec):
                abundances[0,j,k] = cell.abundances[species[k]]

    else:
        
        oldoutput = np.load(f'{outputdirr}/{outputfile}')
        # get things that don't change
        times    = oldoutput['times']
        species  = oldoutput['species']
        heights  = oldoutput['heights']
        gas_dens = oldoutput['density']
        # update the things that do
        extinctions = oldoutput['extinctions']
        abundances  = oldoutput['abundances']
        nspec = len(species)
        for j in range(nzs):
            cell = col.cells[j]
            extinctions[nout,j] = cell.av
            for k in range(nspec):
                abundances[nout,j,k] = cell.abundances[species[k]]

    np.savez(f'{outputdirr}/{outputfile}',times=times,species=species,
        abundances=abundances,heights=heights,density=gas_dens,
        extinctions=extinctions)
    if delete_subdir:
        print('removing sub directories')
        print(f'rm -rf {outputdirr}/z*')
        subprocess.run(f'rm -rf {outputdirr}/z*', shell=True)


def get_nh_and_av(source_file: str) -> tuple[float, float]:
    """Read the number denisty of hydrogen and visual extinction from
    the astrochem source file

    Args:
        source_file (str): name of astrochem source file

    Returns:
        tuple[float, float]: number density of hydrogen and visual extinction
    """
    with open(source_file,'r') as f:
        f.readline()  # skip header
        params = f.readline().split()
    nh = float(params[2])
    av = float(params[1])
    return nh,av

def get_specs_array(outputdirr: str) -> np.ndarray:
    """Returns an array of species used from astrochem output. Reads in
    from outputdirr/z00/astrochem_output.h5

    Args:
        outputdirr (str): name of output directory

    Returns:
        np.ndarray: array of species names
    """

    _,abundict = get_abundict(outputdirr+f'/z00/astrochem_output.h5')
    specs = np.empty(len(list(abundict.keys())),dtype='U25')
    for k,spec in enumerate(list(abundict.keys())):
        specs[k] = spec
    return specs

def read_touts_from_file(fname: str) -> np.ndarray:
    """Helper to read in output times saved in times.out

    Args:
        fname (str): name of file

    Returns:
        np.ndarray: array of output times
    """
    touts = []
    with open (fname,'r') as f:
        f.readline() # header
        for line in f:
            t,time = line.split('\t')
            touts.append(float(time))
    return np.array(touts)

def gather_all_abundances(
        col: Column, outputdirr: str, outputfile: str, 
        tf: float, delete_subdir: bool) -> None:
    """(DEPRECATED) reads in all of the astrochem outputs from a column
    and puts all the output params into a single .npz file

    Args:
        col (Column): CANDY column
        outputdirr (str): output directory
        outputfile (str): output file
        tf (float): final time of simulation
        delete_subdir (bool): option to delete subdirectories after
        gathering abundances
    """
    specs = get_specs_array(outputdirr)
    nspec = len(specs)
    touts = read_touts_from_file(f'{outputdirr}/times.out')
    touts = touts[touts<=tf]
    nts = len(touts)
    nzs = col.ncells

    gas_dens = np.zeros(nzs)
    extinctions = np.zeros((nts,nzs))
    abundances = np.zeros((nts,nzs,nspec))
    heights = np.zeros(nzs)
    for j in range(nzs):
        dirr = f'{outputdirr}/z{j:0>2}'
        for t,time in enumerate(touts):
            nh,av = get_nh_and_av(dirr+f'/source_t{t}.mdl')
            extinctions[t,j] = av
            if t==0:
                gas_dens[j] = nh
                heights[j] = col.cells[j].z
            final_dict = get_final_abuns(dirr+f'/astrochem_output_t{t}.h5')
            for k,spec in enumerate(list(final_dict.keys())):
                abundances[t,j,k] = final_dict[spec]
    np.savez(f'{outputdirr}/{outputfile}',times=touts,species=specs,
        abundances=abundances,heights=heights,density=gas_dens,
        extinctions=extinctions)
    if delete_subdir:
        print('removing sub directories')
        print(f'rm -rf {outputdirr}/z*')
        subprocess.run(f'rm -rf {outputdirr}/z*', shell=True)

def gather_static_abuns(
        col: Column, outputdirr: str, outputfile: str, 
        tf: float, delete_subdir: bool):
    """Gather outputs from column into single .npz file for a static
    run of CANDY

    Args:
        col (Column): CANDY column
        outputdirr (str): output directory
        outputfile (str): output filename
        tf (float): final time from the simulation
        delete_subdir (bool): option to delete subdirectories after
        creating .npz file
    """
    specs = get_specs_array(outputdirr)
    nspec = len(specs)
    times,_ = get_abundict(f'{outputdirr}/z00/astrochem_output.h5')
    nts = len(times)
    nzs = col.ncells

    gas_dens = np.zeros(nzs)
    heights = np.zeros(nzs)
    extinctions = np.zeros((nts,nzs))
    abundances = np.zeros((nts,nzs,nspec))
    for j in range(nzs):
        dirr = f'{outputdirr}/z{j:0>2}'
        nh,av = get_nh_and_av(dirr+'/source.mdl')
        extinctions[:,j] = np.ones_like(times)*av
        gas_dens[j] = nh
        heights[j] = col.cells[j].z
        ts,abundict = get_abundict(dirr+'/astrochem_output.h5')
        for k,spec in enumerate(list(abundict.keys())):
            abundances[:,j,k] = abundict[spec]

    np.savez(f'{outputdirr}/{outputfile}',times=times,species=specs,
        abundances=abundances,heights=heights,density=gas_dens,
        extinctions=extinctions)
    if delete_subdir:
        print('removing sub directories')
        print(f'rm -rf {outputdirr}/z*')
        subprocess.run(f'rm -rf {outputdirr}/z*', shell=True)

def readh5file(filename: str='astrochem_output.h5') -> tuple[h5.File, h5.Dataset, h5.Dataset, h5.Dataset]:
    """Helper function to read astrochem output and return relevent information

    Args:
        filename (str, optional): name of .h5 file to read. Defaults to 'astrochem_output.h5'.

    Returns:
        tuple[h5.File, dict, dict, dict]: tuple of outputs
            opened h5.File object (be sure to close this later!)

            list of species from astrochem

            array of abundances for each species over time

            list of times
    """

    f = h5.File(filename,'r')
    spec = f['Species']
    abun = f['Abundances']
    time = f['TimeSteps']
    return f,spec,abun,time

def get_abundict(filename='astrochem_output.h5',specs = 'all') -> tuple[np.ndarray, dict]:
    """Create a dictionary of species abundances from an astrochem output

    Args:
        filename (str, optional): astrochem output h5 file. Defaults to 'astrochem_output.h5'.
        specs (list[str] | str, optional): list of species to include in the output dictionary. Defaults to 'all'.

    Returns:
        tuple[np.ndarray, dict]: array of times and the dictionary of species and abundances
    """
    f,spec,abun,time = readh5file(filename)

    d = {}
    time = np.array(time)
    for i in range(spec.size):
        specie = spec[i].decode('utf-8')
        if specie in specs or specs == 'all':
            abund = abun[0,:,i]
            d[specie] = abund
    f.close()
    return time,d

def get_final_abuns(f_name: str,specs: list[str]|str='all') -> dict:
    """Create a dictionary of the species abundances at the last time of the astrochem output

    Args:
        f_name (str): name of astrochem .h5 file
        specs (list[str] | str, optional): list of species to include. Defaults to 'all'.

    Returns:
        dict: dictionary of species and abundances
    """
    final_abuns = {}
    times,d = get_abundict(f_name,specs=specs)
    for spec in d:
        final_abuns[spec]=d[spec][-1]
    return final_abuns

def save_pebbles(
        col: Column, pebcomp: dict, f_pebout: str, time: float) -> None:
    """Save the pebble compistion to text file

    Args:
        col (Column): CANDY column
        pebcomp (dict): dictionary of pebble compositions
        f_pebout (str): name of pebble text file
        time (float): the current time
    """
    with open(f_pebout, 'a+') as f:
        for spec in list(pebcomp.keys()):
            peb_comp_norm = pebcomp[spec]/col.cells[0].NH
            f.write(f'{time}\t{spec}\t{peb_comp_norm:.10e}\n')


if __name__ == '__main__':
    import argparse
    from pprint import pprint
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',default='example.ini')
    args = parser.parse_args()
    model,phys,abuns = read_infile(args.infile)
    pprint(model)
    pprint(phys)
    pprint(abuns)