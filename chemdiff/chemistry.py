from asyncio import subprocess
import subprocess
import multiprocessing as mp
import os

from .candyio import get_final_abuns
from .column import Column

def do_chemistry(
        col: Column, chemtime: float, f_chm: str, outdirr: str,
        abs_err=1.e-20, rel_err=1.e-10) -> None:
    """ Setup chem_helper() function to call astrochem in parallel
    using python multiprocessing library.

    Parameters
    ----------
    col
        Column to do chemistry on
    chemtime
        Time (in years) over which to do chemistry
    f_chm
        chm file to use for chemistry
    outdirr
        output directory
    abs_err, rel_err
        absolute and relative errors for chemistry integration
    """
    args = [(col, j, f_chm, chemtime, abs_err, rel_err, outdirr) 
             for j in range(col.ncells)]
    with mp.Pool() as pool:
        solvedcells = pool.map(chem_helper,args)
    for j in range(col.ncells):
        col.cells[j].update_abundances(solvedcells[j])

def chem_helper(args: tuple) -> dict:
    """ Helper function to parallelize chemistry calculation. Calls
    Astrochem on a given cell

    Parameters
    ----------
    args
        tuple of arguments from do_chemistry() function
    
    Returns
    -------
    dict
        update abundance dictionary after chemistry
    """
    cwd = os.getcwd()
    col, j, f_chm, chemtime, abs_err, rel_err, outdirr = args
    dirr = f'{cwd}/{outdirr}/z{j:0>2}'
    cell = col.cells[j]
    cell.write_chem_inputs(chemtime, abs_err, rel_err, f_net=f_chm, 
        f_input=f'{dirr}/input.ini', f_source=f'{dirr}/source.mdl')
    # print('working on cell ',j)
    subprocess.run(['astrochem','-q','input.ini'],cwd=dirr)
    # print('done with cell ',j)
    d = get_final_abuns(f'{dirr}/astrochem_output.h5','all')
    return d

def make_chmfile(abuns: dict) -> str:
    """Make a chm file with each species doing nothing so that their
    abundances are still returned in the output file.

    Args:
        abuns (dict): dictionary where the keys are species to be included in the do-nothing network

    Returns:
        str: name of file created
    """
    # template
    # grain -> grain    0.00e+00  0.00e+00  0.00e+00  2  6830
    reacno = 1
    chmname = 'nochem.chm'
    with open(chmname,'w+') as f:
        for spec in abuns:
            line = f'{spec} -> {spec}    0.00e+00  0.00e+00  0.00e+00  2  {reacno}\n'
            f.write(line)
            reacno+=1
        if 'grain' not in abuns:
            line = f'grain -> grain    0.00e+00  0.00e+00  0.00e+00  2  {reacno}\n'
            f.write(line)
            reacno+=1
    return chmname
