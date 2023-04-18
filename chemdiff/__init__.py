import numpy as np
import h5py as h5
import subprocess

__version__ = "1.1.0"
__author__ = "Eric Van Clepper"

from .column import Column
from .cell import Cell
from . import candyio
from . import constants
from . import disk
from . import chemistry
from . import diffusion
from . import utils

def create_column(r: float, alpha=1e-3, ncells=50) -> Column:
    """Create a Column object.

    Parameters
    ----------
    r
        radial location of column in cm
    alpha
        turbulent alpha parameter
    ncells
        number of cells in the column
    
    Returns
    -------
    Column
    """
    return Column(r,alpha,ncells)

def create_cell(r: float, z: float, **kwargs) -> Cell:
    """ Create a cell object
    
    Parameters
    ----------
    r
        radial location in cm
    z
        elevation of cell in cm
    kwargs
        key word arguments to pass onto cell.Cell object

    Returns
    -------
    Cell
    """
    return Cell(r,z,**kwargs)