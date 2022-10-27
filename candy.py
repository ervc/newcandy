import os, subprocess

import chemdiff as cd
import chemdiff.constants as const

def main(infile):
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

    ### main loop
    if not STATIC:
        nout = 0
        time = model_dict['ti']
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
                        phys_dict['growth_height'])
                subtime += diffdt

            col.update_column_densities(phys_dict['opacity'])

            time += chemdt
            if time >= tout:
                print(f'time: {time}')
                cd.candyio.save_outputs(col,nout,outputdirr)
                if GROWTH:
                    cd.candyio.save_pebbles(col,pebcomp,
                        f'{outputdirr}/{f_pebout}',time)
                nout+=1

        ### gather abundances
        cd.candyio.gather_all_abundances(col,outputdirr,outfile,
            model_dict['tf'],True)
    else:
        print('static run')
        cd.chemistry.do_chemistry(col,chemtime,f_chm,outputdirr)
        cd.candyio.gather_static_abuns(col,outputdirr,outfile,
            model_dict['tf'],True)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    args = parser.parse_args()
    main(args.infile)