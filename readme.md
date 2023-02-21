# Candy

Chemistry ANd DYnamics in protoplanetary disks

- [Candy](#candy)
  - [About](#about)
  - [Installation and Requirements](#installation-and-requirements)
  - [Quickstart](#quickstart)
    - [Static runs](#static-runs)
    - [Diffusive runs](#diffusive-runs)
    - [Growth runs](#growth-runs)
  - [Input parameters](#input-parameters)
    - [Model](#model)
    - [Phys](#phys)
    - [Abundances](#abundances)
  - [Output files](#output-files)
    - [outfile](#outfile)
    - [Pebble composition](#pebble-composition)
  - [Astrochem Changes](#astrochem-changes)
  - [chemdiff](#chemdiff)

## About

CANDY (Chemistry ANd DYnamics) is a protoplanetary disk model for determining the time-dependent
chemistry and a 1D vertical slice of a dynamic disk. CANDY includes
chemistry, vertical diffusion, pebble growth, and ice sequestration.
Details of the model are given in [Van Clepper et al. 2022](https://iopscience.iop.org/article/10.3847/1538-4357/ac511b).

CANDY contains two submodules, a modified astrochem ([Maret & Bergin
2015](http://ascl.net/1507.010)) directory and the chemdiff directory, where the python wrapper is
contained. CANDY works by iteratively switching between chemistry and dynamics for a vertical Column at a given location in a protoplanetary disk. The column is split into adjacent cells, within which the chemical and physical properties are assumed to be constant. Each is allowed to evolve chemically independently for some period of time before gas and ice species diffuse vertically through the disk between adjacent cells. Meanwhile, small grains grow into larger pebbles, removing ice species from the active chemistry.


## Installation and Requirements

Astrochem dependencies are given in the astrochem wiki, while chemdiff dependencies include openmpi and mpi4py. Recommended installation is using conda to create a virtual environment with all dependencies. A new environment can be created with all needed dependencies with

	conda create --name newcandy -c conda-forge sundials=5.7.0 python cython numpy matplotlib h5py autoconf automake libtool mpi4py openmpi

Active the enviornment with

	conda activate newcandy

or

	source activate newcandy

To compile Astrochem run the following commands *from the Astrochem subdirectory*

	./bootstrap
	./configure CPPFLAGS="-I$CONDA_PREFIX/include" LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib -L/$CONDA_PREFIX/lib"  --prefix=$CONDA_PREFIX
	make
	make install

Because candy.py is run from the parent directory, there is no need to adjust your python path!

## Quickstart

To run candy, simply run 

    python candy.py examples/growth_example/growth_example.ini

from the command line, where growth_example.ini is the input file you wish to use.

### Static runs

In a static candy run, there is no diffusion nor is there any grain growth. 

To set up a static run of candy, set

    difftime = 0

in the model parameters of the input file. Additionally, the `chemtime` and `tf` parameters must be set equal to one another. This is because in a static run, chemistry is called only once per cell. To avoid ambiguity in how long to run chemistry for, these parameters should be the same.

### Diffusive runs

In a diffusive run, there is diffusion of material throughout the column, but grains do not grow into pebbles. This can be done by setting the `difftime` and `chemtime` model parameters, but the physical parameter `growth_height` should be set to

    growth_height = 0

This sets the height below which pebbles will grow to 0, shutting off any pebble growth.

*Note:* for diffusive and growth runs, diffusion is the time limiting step. Therefore we must have

    difftime <= chemtime <= tf

### Growth runs

A growth run contains diffusion and pebble growth (and is the default setup for CANDY). The `growth_height` parameter sets the height (in disk scaleheight) below which pebbles will grow. Any[^1] `growth_height` > 0 will result in pebble growth. The resulting pebble composition is saved in the `pebfile` (default: `pebcomp.out`) in the `outputdir` directory.

Pebble compositions are stored normalized to the column density of hydrogen, with units of cm$^{-2}$.

[^1]: The `growth_height` must be larger than the z-value of the lowest cell in the column. The column is 5 scaleheights tall, split between `ncells` cells. Cell z-values are stored for the center of each cell, so `growth_height` must be above $0.5*5/$`ncells`.

## Input parameters

Input files for CANDY are split into three sections: model inputs `[model]` which controls solver parameters, physical inputs `[phys]` which sets the relevent physical values for chemistry, and the molecular abundances `[abundances]`, which sets the intial abundances of each chemical species in the cells.

Any parameters that are not explicitly set in the input file will use the default values.

*Note:* Any relative paths to files are always from the parent directory (where `candy.py` is located). For example, if you had an input file in an input directory (say `newcandy/inputs/candyinput.ini`), the path to the chmfile will be relative to the `newcandy` directory. Output files will automatically be placed into the output directory `outputdir`.

### Model

- chmfile
  - Chemical network file to be used by astrochem. *Default:* astrochem/networks/umist12_x.chm
- pebfile
  - name of file where pebble compositions will be stored. *Default:* pebcomp.out
- ti
  - Start time of candy integration in years. *Default:* 0
- tf
  - End time of the candy integration in years. *Default:* 1.0e6
- ncells
  - Number of cells to split the column into. *Default:* 50
- chemtime
  - Maximum time in years to run chemistry in each cell before diffusion. *Default:* 500
- difftime
  - Maximum timestep for diffusion between cells and pebble growth. *Default:* 100
- outputdir
  - Output directory where output files will be stored (trailing `/` is not needed, this is always interpreted as a directory).  *Default:* r00
- outfile
  - Name of npz file where all outputs will be saved. *Default:* output_abuns.npz

### Phys

- r
  - radial location of integration. *Default:* 30
- r_units
  - units of radius, options are 'au' or 'cm'. *Default:* au
- alpha
  - Alpha viscosity parameter of the disk. *Default:* 1.0e-3
- chi
  - external UV radiation field strength in Draine Units. *Default:* 50
- cosmic
  - Cosmic ionization rate in s$^{-1}$. *Default:* 1.3e-17
- grain_size
  - size of solid grains in microns. *Default:* 0.1
- dg0
  - Initial dust-to-gas mass ratio. *Default:* 0.01
- opacity
  - The UV opacity of the dust in cm$^2$ g$^{-1}$. *Default:* 1.0e5
- growth_timescale_factor
  - Sets the timescale of growth for small grains into large pebbles. The growth timescale is defined as $\tau_\mathrm{grow} = a/\Omega\epsilon$, where a is the growth_timescale_factor. *Default:* 1
- growth_height
  - Scaleheights below which pebbles will form. *Default:* 1
- growth_delay_time
  - Time in years to wait before beginning pebble growth. Grains and ice will grow into pebbles when time >= growth_delay_time. *Default:* 0

### Abundances

The initial abundances of each molecular species. Abundances are assumed to be constant throughout the column intially, default abundances are the same as listed in table 1 of Van Clepper et al. 2022, reproduced below:

| Molecule | Abundance |
|---------:|:----------|
H$_2$ | 5.00(-1)
He | 9.75(-2)	
NH$_3$ | 1.45(-6)
H$_2$O | 1.18(-4)	
CO | 6.00(-5) 
N$_2$ | 2.00(-5)	
CH$_4$ | 2.00(-6) 
CH$_3$OH | 1.00(-6)	
H$_2$S | 1.91(-8) 
CO$_2$ | 5.00(-5)	
HCN | 3.50(-7) 
Grains | 2.20(-12)

## Output files

As candy runs, it will create subdirectories named `z00/`, `z01/`, etc... for each cell in the Column. Within each of these folders will be `astrochem_output` HDF5 files for each output time from the simulation. The default behavior is for these many files to be deleted at the end of the CANDY run, and all of the outputs are gathered into one `.npz` file. CANDY also creates two plain text files: the pebble composition (default `pebcomp.out`) and output times, `times.out`.

### outfile

The main output file, outfile in the [model parameters](#model), is saved as a zipped, uncompressed `.npz` file. This file can be opened in Python, by importing the `numpy` package and using the `np.load()` function (see documentation [here](https://numpy.org/doc/stable/reference/generated/numpy.load.html#numpy.load)). For the default outfile name of 'output_abuns', the file could be loaded as follows:

```python
import numpy as np
cdout = np.load('output_abuns.npz')
```

This returns an NpzFile object, containing the saved numpy arrays that are output from CANDY. The available array names can be found by accessing the `files` property

```python
cdout.files
# out: ['times', 'species', 'abundances', 'heights', 'density', 'extinctions']
```

These arrays can be indexed similar to a python dictionary, using square brackets. For example, 'times' contains a list of output times.

```python
times = cdout['times']
times
# out: array([  100.,   200.,   500.,  1000.,  2000.,  5000., 10000.])
```

The different arrays have different shapes dependant on which dimension they vary. For an output with `nt` times, `nz` cells, and `nspec` species, the outputs are as follows:

- times: *1D array, shape (nt,)*
  - Output times in years.
- species: *1D array, shape (nspec,)*
  - List of all species considered in the network.
- abundances: *3D array, shape (nt, nz, nspec)*
  - Abundances relative to hydrogen of all species, at all heights and all times.
- heights: *1D array, shape (nz,)*
  - The elevation above the midplane at the center of each cell in cm.
- density: *1D array, shape (nz,)*
  - The number density of hydrogen in a given cell in cm$^{-3}$.
- extinctions: *2D array, shape (nt, nz)*
  - The visual extinction, $A_v$ of a given cell at a given time. Convert to the UV optical depth using the relation $A_v = \tau_{UV}/3.02$
  

See below for a (quick) example of how to read abundances, this is also demonstrated in an example in the `plot_results.py` plotting script.

```python
import numpy as np
# load data
cdout = np.load('examples/growth_example/output_abuns.npz')
times = cdout['times']
specs = cdout['species']
abuns = cdout['abundances']
elevs = cdout['heights']
# Say we want to find the midplane abundance of H2O ice after 1000 years of evolution
# find the correct indices using np.where() function
t = np.where(times == 1000)[0][0]       # time = 1000 years
j = 0                                   # midplane, cell index = 0
k = np.where(specs == 'grainH2O')[0][0] # species index of water ice
abuns[t,j,k]
# out: 1.6693156020606896e-05
```

### Pebble composition

The pebble composition file is saved as plain text, default name `pebcomp.out`. This file keeps track of the ice species that are sequestered away from the active chemistry on the mantle of growing pebbles. The abundances of each ice species are normalized to the column density of hydrogen, rather than the number density of hydrogen. This is because the entire column will only have *one* reservoir of pebbles. The routine for removing ice is in the `grow_grains()` function in the `chemdiff/diffusion.py` module, and this is written to the pebcomp.out file in the `save_pebbles()` function in `chemdiff/candyio.py`.

The pebble composition file has three columns, the first is the time of the output in years, second the name of the ice species, and the third the normalized abundance.

## Astrochem Changes

The astrochem folder is copied from the [astrochem](https://github.com/smaret/astrochem) github page, complete with documentation and installation instructions. Changes have been made within the `src` folder, and Astrochem can be installed as usual. If you already have astrochem installed on your machine, you can simply copy the `src` folder from here into your astrochem directory (replacing the default astrochem/src folder), then reinstall as usual.

Chemical network files (`.chm`) should follow the same format as astrochem but with added chemical reactions for:

|description | reacno|
|------------:|:-------|
| HD formation on grains | 100 |
| Shielded dissociation of H2 | 14 |
| Shielded dissociation of HD | 114 |
| Shielded photo-ionization of C | 15 |
| Shielded dissociation of C16O | 16 |
| Shielded dissociation of C17O | 17 |
| Shielded dissociation of C18O | 18 |
| Hydrogenation | 25 |
| Cosmic-ray desorption of CO | 42 |
| Cosmic-ray reactions involving He | 43 |
| secondary xray ionization of ... H | 60 |
| ... H2 | 61 |
| ... other molecules | 62 |
| creation of excited H2 | 90 |
| de-excitation of H2 | 91 | 
| reactions with excited H2 | 92 |

Some reaction rate calculations are also adjusted to be more similar to the DALI calculations (Burderer et al 2012). A full description of the handling of the chemical calculations is in the appendix of Van Clepper et al. 2022.

## chemdiff

The python wrapper for for running astroCHEM with DIFFusion. This is entirely written in python, and successively calls `astrochem` in parallel at different vertical locations ranging from disk midplane to 5 scale heights. Input parameters are readin from custom `.ini` files -- examples have been provided for growth, static, and diffusion setups.

The chemdiff directory contains all of the python modules. Of particular interest:
- `chemistry.py`
  - Calls astrochem in parallel for each of the cells.
- `diffusion.py`
  - Diffuses species between cells and contains grain growth functions.
- `disk.py`
  - Contains the disk parameters, including surface density and temperature profiles.
- `utils.py`
  - Contains miscellanious, possibly useful, functions for analysis.