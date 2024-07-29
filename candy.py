import os, subprocess

import chemdiff as cd
import chemdiff.constants as const

def main(infile,CONT=-1):
    ''' Main CANDY algorithm, switch off between chemistry and dynamics
    Parameters
    ----------
    infile : str
        name of input file to read parameters from
    CONT : int, default=-1
        if greater than zero, continue run from the CONT output. 
        Otherwise, create the cells using abundances from the infile

    '''
    print(f'*** NEWCANDY v{cd.__version__} ***',flush=True)
    ### readin inputs
    model_dict,phys_dict,init_abuns = cd.candyio.read_infile(infile)
    STATIC = False
    GROWTH = True
    if model_dict['difftime'] <= 0:
        STATIC = True
        GROWTH = False
        phys_dict['growth_height'] = 0
    elif phys_dict['growth_height'] <= 0:
        GROWTH = False

    if STATIC:
        assert model_dict['tf'] == model_dict['chemtime'], "For static run (difftime = 0), chemtime must be same as tf"


    ### setup the output files
    cwd = os.getcwd()
    f_pebout = model_dict['pebfile']
    if model_dict['chmfile'] == 'None':
        f_chm = cwd+'/'+cd.chemistry.make_chmfile(init_abuns)
    else:
        f_chm = cwd+'/'+model_dict['chmfile']
    outputdirr = model_dict['outputdir']
    outfile = model_dict['outfile']

    difftime = float(model_dict['difftime'])
    chemtime = float(model_dict['chemtime'])
    touts = model_dict['touts']

    if not os.path.exists(outputdirr):
        subprocess.run(['mkdir',outputdirr])

    if GROWTH:
        pebcomp = {}
        with open(f'{outputdirr}/{f_pebout}','w+') as f:
            f.write('Time [yr]\tSpecies\tPebble Abundance [cm-2]\n')

    if not STATIC:
        with open(f'{outputdirr}/times.out','w+') as f:
            f.write('output number\toutput time [yr]\n')
            for t,tout in enumerate(touts):
                f.write(f'{t}\t{tout}\n')

    ### initialize the column
    r = phys_dict['r']
    if phys_dict['r_units'].lower() == 'au':
        r*=const.AU
    elif phys_dict['r_units'].lower() != 'cm':
        raise NameError("r_units must be au or cm, default is au")

    alpha = float(phys_dict['alpha'])
    nzs = model_dict['ncells']

    col = cd.create_column(r, alpha, nzs)

    if (not STATIC) and (CONT<0):
        if 'grain' not in init_abuns:
            print('grain abundance is not given, creating grain abundance from dg0')
            grabun = cd.utils.dg2grain_abun(phys_dict['dg0'],phys_dict['grain_size'])
            print(f'grain abundance = {grabun:.4e}',flush=True)
            init_abuns['grain'] = grabun
        else:
            print('grain abundance is given, please make sure grain species is in the chemical network somewhere!')
            dg_from_abun = cd.utils.grain_abun2dg(init_abuns['grain'])
            print(f'dust-to-gas ratio from grain abundance = {dg_from_abun:.2e}')
            print(f'dust-to-gas ratio from input file = {phys_dict["dg0"]:.2e}')
            print('Are these consistent???',flush=True)

    if CONT < 0:
        ### create the cells
        ### integrate column densities while we do this
        ### calculating column density from top down
        NCO = 0.
        NH2 = 0.
        NHD = 0.
        NH = 0.
        tau = 0.
        xray = 0.
        for j in reversed(range(nzs)):
            dirr = f'{outputdirr}/z{j:0>2}'
            if not os.path.exists(dirr):
                subprocess.run(['mkdir',dirr])
            # move the data_coselfs_iso.dat file into directory where
            # astrochem will be run
            subprocess.run(['cp','data_coselfs_iso.dat',
                f'{dirr}/data_coselfs_iso.dat'])
            z = col.dz*(j+0.5)
            rho = cd.disk.get_density(r,z)
            temp = cd.disk.get_temperature(r,z)
            nh = 2*rho/const.MBAR
            nco = 0.
            if 'CO' in init_abuns:
                nco = nh*init_abuns['CO']
            nh2 = 0.
            if 'H2' in init_abuns:
                nh2 = nh*init_abuns['H2']
            nhd = 0
            if 'HD' in init_abuns:
                nhd = init_abuns['HD']*nh
            NCO += nco*col.dz
            NH2 += nh2*col.dz
            NHD += nhd*col.dz
            NH += nh*col.dz
            # optical depth
            tau += rho*phys_dict['opacity']*col.dz*phys_dict['dg0']

            col.cells[j] = cd.create_cell(
                r,z,chi=phys_dict['chi'],
                cosmic=phys_dict['cosmic'],grain_size=phys_dict['grain_size'],
                dust_gas_ratio=phys_dict['dg0'],av=tau/3.02,rho=rho,
                Tgas=temp,Tdust=temp,xray=xray,NCO=NCO,NH2=NH2,NHD=NHD,NH=NH,
                abundances=dict(init_abuns))
    else:
        ### Read in cells from CONT output
        NH = 0.
        for j in reversed(range(nzs)):
            dirr = f'{outputdirr}/z{j:0>2}'
            # get some params from cdinput
            z = col.dz*(j+0.5)
            # read in files
            fastro = dirr+f'/astrochem_output_t{CONT}.h5'
            fsource= dirr+f'/source_t{CONT}.mdl'
            finput = dirr+f'/input_t{CONT}.ini'
            # get abundances from last astro out
            abundict = cd.candyio.get_final_abuns(fastro)
            # get params from source file
            with open(fsource,'r') as f:
                # header
                f.readline()
                (n,av,nh,tgas,tdust,NCO,NH2,NHD,
                    xray,rau,zau) = list(map(float,f.readline().split()))
            rho = nh/2*const.MBAR
            NH += nh*col.dz
            dg = phys_dict['dg0']
            with open(finput,'r') as f:
                for line in f:
                    if line.startswith('grain_gas_mass_ratio'):
                        dg = float(line.split()[-1])
            col.cells[j] = cd.create_cell(
                r,z,chi=phys_dict['chi'],
                cosmic=phys_dict['cosmic'],grain_size=phys_dict['grain_size'],
                dust_gas_ratio=dg,av=av,rho=rho,
                Tgas=tgas,Tdust=tdust,xray=xray,NCO=NCO,NH2=NH2,NHD=NHD,NH=NH,
                abundances=dict(abundict))


    ### main loop
    if not STATIC:
        if CONT<0:
            print('Starting new run')
            nout = 0
            time = model_dict['ti']
        else:
            nout = CONT+1
            time = touts[CONT]
            print(f'Continuing from output {CONT}, current time = {time}')
        tf = model_dict['tf']
        while time < tf:
            tout = touts[nout] # next time to output results for
            print('next tout: ',tout)
            chemdt = min(chemtime,tout-time)

            #### DO CHEMISTRY ####
            sols = cd.chemistry.do_chemistry(col, chemdt, f_chm, outputdirr)

            #### DO DIFFUSION AND GROWTH ####
            subtime = 0
            while (subtime < chemdt):
                diffdt = min(difftime,chemdt-subtime)

                cd.diffusion.do_diffusion(col,diffdt)
                if GROWTH and time >= phys_dict['growth_delay_time']:
                    pebcomp = cd.diffusion.grow_grains(col,diffdt,pebcomp,
                        phys_dict['growth_timescale_factor'],
                        phys_dict['growth_height'],
                        outputdirr)
                subtime += diffdt

            col.update_column_densities(phys_dict['opacity'])

            ### step time and save values ###
            time += chemdt
            if time >= tout:
                print(f'time: {time}')
                dsd = False
                if time >= tf: dsd = True
                cd.candyio.new_save_outputs(col,nout,time,
                    outputdirr,outfile,delete_subdir=dsd)
                if GROWTH:
                    cd.candyio.save_pebbles(col,pebcomp,
                        f'{outputdirr}/{f_pebout}',time)
                nout+=1

    else:
        print('static run')
        cd.chemistry.do_chemistry(col,chemtime,f_chm,outputdirr)
        cd.candyio.gather_static_abuns(col,outputdirr,outfile,
           model_dict['tf'],True)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("-C","--Continue",default=-1,type=int)
    args = parser.parse_args()
    main(args.infile,args.Continue)