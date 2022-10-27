/*
   rates.c - Compute the reaction rates.

   Copyright (c) 2006-2019 Sebastien Maret

   This file is part of Astrochem.

   Astrochem is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   Astrochem is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Astrochem.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrochem.h"

double
getcoss (double NCO, double NH2, int isotope)
{
  char cdummy[800];
  int nNCO, nNH2, i, j;
  double theta = 1.;
  FILE *fp;

  if (isotope != 16 && isotope != 17 && isotope != 18) {
    fprintf (stderr, "Unknown isotope of Carbon %i, known isotopes of carbon are 16, 17, and 18",isotope);
    exit(1);
  }
  // data_coselfs_ios.data should be in project directory, but
  // astrochem is run from project/r00/zXX, so look up two directories!
  fp=fopen("../../data_coselfs_iso.dat","r");
  if(fp==NULL){ printf("STOP: Cannot read data_coselfs_iso.dat\n"); exit(1); } 

  // fprintf(stdout, "Height is %f AU\n= %f cm\n", zau,zcm); 

  for(i=0;i<5;i++) fgets(cdummy, 800, fp);
  fgets(cdummy, 24, fp); fscanf(fp,"%i\n", &nNCO);
  fgets(cdummy, 24, fp); fscanf(fp,"%i\n", &nNH2);

  if(nNCO!=47 || nNH2!=42){
    printf("STOP: Format problem in data_coselfs_iso.dat\n"); exit(1);  
  }

  // create arrays for the list of NCO, NH2, and arrays for C16O, C17O, C18O factors
  double NCOs[nNCO], NH2s[nNH2];
  double c16oss_factors[nNCO][nNH2], c17oss_factors[nNCO][nNH2], c18oss_factors[nNCO][nNH2];

  // store NCO, NH2 values from file as log 10 values
  fgets(cdummy, 800, fp); 
  for(i=0;i<nNCO;i++){
    fscanf(fp,"%lf\n",&(NCOs[i]));
    NCOs[i]=log10(NCOs[i]);
  }
  
  fgets(cdummy, 800, fp); 
  for(i=0;i<nNH2;i++){
   fscanf(fp,"%lf\n",&(NH2s[i]));    
    NH2s[i]=log10(NH2s[i]);   
  }

  // read in self shielding functions for isotopes of oxygen
  // note: format of data_coselfs_iso.dat has been changed from dali
  //       to only include O-isotope data
  fgets(cdummy, 800, fp);
  for(j=0;j<nNH2;j++){ 
    for(i=0;i<nNCO;i++){ 
      fscanf(fp,"%lf\n",&(c16oss_factors[i][j]));    
      c16oss_factors[i][j]=log10(c16oss_factors[i][j]);
    }
  }

  fgets(cdummy, 800, fp);
  for(j=0;j<nNH2;j++){ 
    for(i=0;i<nNCO;i++){ 
      fscanf(fp,"%lf\n",&(c18oss_factors[i][j]));    
      c18oss_factors[i][j]=log10(c18oss_factors[i][j]);
    }
  }

  fgets(cdummy, 800, fp);
  for(j=0;j<nNH2;j++){ 
    for(i=0;i<nNCO;i++){ 
      fscanf(fp,"%lf\n",&(c17oss_factors[i][j]));    
      c17oss_factors[i][j]=log10(c17oss_factors[i][j]);
    }
  }

  fclose(fp);

  double collCO = log10(fmax(1.0,NCO)); // if there is no CO or H2 then set to 1 to avoid issues with log
  double collH2 = log10(fmax(1.0,NH2));

  // find the closest (without going over) column densities from given list
  int icollCO, jcollH2;
  for (i=0;i<nNCO;i++) {
    if (collCO>=NCOs[i]) {
      icollCO = i;
    }
  }
  for (j=0;j<nNH2;j++) {
    if (collH2>=NH2s[j]) {
      jcollH2 = j;
    }
  }

  // fprintf(stdout, "icollCO = %d, jcollH2 = %d \n", icollCO,jcollH2);
  // fprintf(stdout, "NCOs[i] = %.1f, NH2s[j] = %.1f\n", NCOs[icollCO],NH2s[jcollH2]);

  double R0,R1,factor;
  // finally, return the self shielding function value
  switch (isotope) {
    case 16:
      // interpolate!
      if (icollCO == nNCO-1){ // max CO, no interpolation in CO direction
        if (jcollH2 == nNH2-1){ // max CO and H2, no interp
          factor = c16oss_factors[icollCO][jcollH2];
          theta = pow(10.0,factor);
          // fprintf(stdout, "No interpolation, factor = %.4f, theta = %.2e\n", factor,theta);
        } else { // max CO only, interpolate in H2 direction
          factor = c16oss_factors[icollCO][jcollH2] * (NH2s[jcollH2+1]-collH2)/(NH2s[jcollH2+1]-NH2s[jcollH2])
               + c16oss_factors[icollCO][jcollH2+1] * (collH2-NH2s[jcollH2])  /(NH2s[jcollH2+1]-NH2s[jcollH2]);
          // fprintf(stdout, "1D H2 interpolation, factor = %.4f\n", factor);
          theta = pow(10.0,factor);
        } 
      } else if (jcollH2 == nNH2-1){ // max H2 but not CO, interp in 1D
        factor = c16oss_factors[icollCO][jcollH2] * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
             + c16oss_factors[icollCO+1][jcollH2] * (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        // fprintf(stdout, "1D CO interpolation, factor = %.4f\n", factor);
        theta = pow(10.0,factor);
      } else {
        R0 = c16oss_factors[icollCO][jcollH2]    * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
          +  c16oss_factors[icollCO+1][jcollH2]  * (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        R1 = c16oss_factors[icollCO][jcollH2+1]  * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
          +  c16oss_factors[icollCO+1][jcollH2+1]* (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        factor = R0 * (NH2s[jcollH2+1]-collH2)/(NH2s[jcollH2+1]-NH2s[jcollH2])
               + R1 * (collH2-NH2s[jcollH2])  /(NH2s[jcollH2+1]-NH2s[jcollH2]);
        // fprintf(stdout, "2D interpolation, factor = %.4f\n", factor);
      
        theta = pow(10.0,factor);
      }
      
      // fprintf(stdout,"Theta c-16-o: %10.3e\n", theta);
      break;
    case 17:
      // interpolate!
      if (icollCO == nNCO-1){
        if (jcollH2 == nNH2-1){
          factor = c17oss_factors[icollCO][jcollH2];
          theta = pow(10.0,factor);
        } else {
          factor = c17oss_factors[icollCO][jcollH2] * (NH2s[jcollH2+1]-collH2)/(NH2s[jcollH2+1]-NH2s[jcollH2])
               + c17oss_factors[icollCO][jcollH2+1] * (collH2-NH2s[jcollH2])  /(NH2s[jcollH2+1]-NH2s[jcollH2]);
          theta = pow(10.0,factor);
        }
      } else if (jcollH2 == nNH2-1){ // max H2 but not CO, interp in 1D
        factor = c17oss_factors[icollCO][jcollH2] * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
             + c17oss_factors[icollCO+1][jcollH2] * (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        theta = pow(10.0,factor);
      } else {
        R0 = c17oss_factors[icollCO][jcollH2]    * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
          +  c17oss_factors[icollCO+1][jcollH2]  * (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        R1 = c17oss_factors[icollCO][jcollH2+1]  * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
          +  c17oss_factors[icollCO+1][jcollH2+1]* (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        factor = R0 * (NH2s[jcollH2+1]-collH2)/(NH2s[jcollH2+1]-NH2s[jcollH2])
               + R1 * (collH2-NH2s[jcollH2])  /(NH2s[jcollH2+1]-NH2s[jcollH2]);
      
        theta = pow(10.0,factor);
      }
      
      // fprintf(stdout,"Theta c-17-o: %e\n", theta);
      break;
    case 18:
      // interpolate!
      if (icollCO == nNCO-1){
        if (jcollH2 == nNH2-1){
          factor = c18oss_factors[icollCO][jcollH2];
          theta = pow(10.0,factor);
        } else {
          factor = c18oss_factors[icollCO][jcollH2] * (NH2s[jcollH2+1]-collH2)/(NH2s[jcollH2+1]-NH2s[jcollH2])
               + c18oss_factors[icollCO][jcollH2+1] * (collH2-NH2s[jcollH2])  /(NH2s[jcollH2+1]-NH2s[jcollH2]);
          theta = pow(10.0,factor);
        }
      } else if (jcollH2 == nNH2-1){ // max H2 but not CO, interp in 1D
        factor = c18oss_factors[icollCO][jcollH2] * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
             + c18oss_factors[icollCO+1][jcollH2] * (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        theta = pow(10.0,factor);
      } else {
        R0 = c18oss_factors[icollCO][jcollH2]    * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
          +  c18oss_factors[icollCO+1][jcollH2]  * (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        R1 = c18oss_factors[icollCO][jcollH2+1]  * (NCOs[icollCO+1]-collCO)/(NCOs[icollCO+1]-NCOs[icollCO])
          +  c18oss_factors[icollCO+1][jcollH2+1]* (collCO-NCOs[icollCO])  /(NCOs[icollCO+1]-NCOs[icollCO]);
        factor = R0 * (NH2s[jcollH2+1]-collH2)/(NH2s[jcollH2+1]-NH2s[jcollH2])
               + R1 * (collH2-NH2s[jcollH2])  /(NH2s[jcollH2+1]-NH2s[jcollH2]);
      
        theta = pow(10.0,factor);
      }
      
      // fprintf(stdout,"Theta c-18-o: %e\n", theta);
      break;
    default:
      theta = 1.;
      fprintf(stderr, "Something went wrong getting the shelf shielding factor\n");
      fprintf(stderr, "check the file %s:%d\n", __FILE__,__LINE__);
      exit(1);
      break;
  }
  return theta;
}


double
rate (double alpha, double beta, double gamm, int reaction_type,
      int reaction_no __attribute__ ((unused)), double nh, double NCO, double NH2, double NHD,
      double av, double tgas, double tdust, double chi, double cosmic,
      double grain_size, double grain_abundance, double ice_abundance, double n_hydro, double n_ice, double xray)
{
  double k;                     /* Reaction rate (cm^-3 s^-1) */
  double theta=1;  /* Self Shielding functions for O isotopes */
  double cosx = cosmic+xray; /*total ionization rate from xrays and CR */
  // set theta to 1 to start, they'll be changed if self shielding is used
  FILE* ff;
  /* Reactions types and rates. The nomencalture is similar to the one
     of the Ohio State University database for astrochemistry, but it
     is extended to include depletion and desorption on/from the grain
     surfaces. */
  double N_des = 2.;
  double dg100 = 100*grain_abundance*GRAIN_MASS_DENSITY_DEFAULT*(4./3.)*M_PI*pow(grain_size,3.)/(CONST_CGSM_MASS_PROTON);
  switch (reaction_type)
    {
    case -1:
      {
        /* Gas-grain interaction (excluding depletion and desorption),
         Electron-grain recombination. */
        k = alpha * pow (tgas / 300, beta) * GAS_DUST_NUMBER_RATIO;
        break;
      }

    case 0:
      {
        /* H2 formation on grains */
        // Cazaux & Tielens 2002
        // k = alpha * pow (tgas / 300, beta) * nh;

        double sqterm, beta_H2, beta_alpha, xi, f_mlps, eta, stick, s_eta;   
        f_mlps=1e-10; // monolayers per second                  
        sqterm=pow(1.0+sqrt((10000.0-200.0)/(600.0-200.0)),2.0); // square root term in Eq. 16 & 17
        beta_H2=3e12*exp(-320.0/tdust); // after Eq. 5
        beta_alpha=0.25*sqterm*exp(-200/tdust); // Eq. 17 
        xi=1.0/(1.0+1.3e13/(2.0*f_mlps)*exp(-1.5*10000.0/tdust)*sqterm); // Eq. 16 (Erratum)
        eta=xi/(1.0+0.005*f_mlps/(2.0*beta_H2)+beta_alpha); // Eq. 15
        
        if(tdust<10.0) eta=1.0; // hack, since Tdust is not realistic for this model
        
        stick=1.0/(1.0+0.04*pow(tgas+tdust,0.5)+2e-3*tgas+8e-6*pow(tgas,2.0));
        s_eta = stick*eta;

        k = s_eta*alpha*pow(tgas,beta)*nh*dg100;
        fprintf(stdout, "k H2 formation %e\n", k);
        break;
      }

    case 100:
      {
        /* HD formation on grains */
        // same as H2 formation, but using H and D abundances in solve.c
        double sqterm, beta_H2, beta_alpha, xi, f_mlps, eta, stick, s_eta;   
        f_mlps=1e-10; // monolayers per second                  
        sqterm=pow(1.0+sqrt((10000.0-200.0)/(600.0-200.0)),2.0); // square root term in Eq. 16 & 17
        beta_H2=3e12*exp(-320.0/tdust); // after Eq. 5
        beta_alpha=0.25*sqterm*exp(-200/tdust); // Eq. 17 
        xi=1.0/(1.0+1.3e13/(2.0*f_mlps)*exp(-1.5*10000.0/tdust)*sqterm); // Eq. 16 (Erratum)
        eta=xi/(1.0+0.005*f_mlps/(2.0*beta_H2)+beta_alpha); // Eq. 15
        
        if(tdust<10.0) eta=1.0; // hack, since Tdust is not realistic for this model
        
        stick=1.0/(1.0+0.04*pow(tgas+tdust,0.5)+2e-3*tgas+8e-6*pow(tgas,2.0));
        s_eta = stick*eta;

        k = s_eta*alpha*pow(tgas,beta)*nh*dg100;
        fprintf(stdout, "k HD formation %e\n", k);
        break;
      }

    case 1:
      {
        /* Cosmic-ray ionization (direct process). Cosmic-ray induced
                   photoreactions (indirect process)   */
        k = alpha * cosmic / COSMIC_DEFAULT;
        break;
      }

    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    {
      /* Ion-molecule reactions, charge exchange reactions (2), Negative
         ion - neutral species reactions (3), Radiative association (4),
         Associative ejection (5), Neutral + neutral -> ion + electron
         (6), Neutral-neutral chemical reactions (7), Neutral-neutral
         radiative association (8), Dissociative recombination (9),
         Radiative recombination (10), Positive ion - negative ion
         recombination (11), Electron attachment (12), Others (14) */
      k = alpha * pow (tgas / 300, beta) * exp (-gamm / tgas);
      break;
    }

    case 13:
    {
      /* Photo-ionization, Photo-dissociation with no shielding*/
      k = alpha * chi * exp (-gamm * av);
      break;
    }

    case 14:
      /* photo-dissociation of H2 with shielding (Draine & Bertoldi eq. 37)*/
      {
        double NH2_5e14, delta_v;
        NH2_5e14 = NH2 / 5.e14;
        delta_v = 0.2; //thermal dispersion in km/s
        theta = 0.965/pow(1.0 + (NH2_5e14/delta_v), 2.0) +
                  0.035/sqrt(1.0 + NH2_5e14) * exp(-8.5e-4*sqrt(1.0 + NH2_5e14));
        k = alpha * chi * exp (-gamm * av) * theta;
        break;
      }
    case 114:
      /* photo-dissociation of HD with shielding (Same as H2 shielding)*/
      {
        double NHD_5e14, delta_v;
        // NHD = 1.e15; // should update from func call!
        NHD_5e14 = NHD/5e14;
        delta_v = 0.2; //thermal dispersion in km/s
        theta=0.965/pow(1.0+(NHD_5e14/delta_v),2.0)+
                0.035/sqrt(1.0+NHD_5e14)*exp(-8.5e-4*sqrt(1.0+NHD_5e14)); 
        theta = 1.;

        k = alpha * chi * exp(-gamm * av) * theta;
        fprintf(stdout, "NHD is %e\n", NHD);
        fprintf(stdout, "k HD -> H+D = %e\n", k);
        break;
      }

    case 15:
    {
      /* Photo ionization of C with shielding factor */
      /* column density of C assumed to be one (no shielding contribution from C, only H2) */
      double NC = 1.;
      theta = exp(-NC*1.1e-17)*
        fmax(0.5,exp(-0.9*pow(tgas,0.27)*pow(NH2/1e22,0.45)));
      k = alpha * chi * exp(-gamm*av) * theta;
      break;
    }

    case 16:
    {
      /* C16O Photo-ionization, Photo-dissociation with shielding */
      theta = getcoss(NCO, NH2, 16);
      k = alpha * chi * exp (-gamm * av) * theta;
      break;
    }

    case 17:
    {
      /* C17O Photo-ionization, Photo-dissociation with shielding */
      theta = getcoss(NCO, NH2, 17);
      k = alpha * chi * exp (-gamm * av) * theta;
      break;
    }

    case 18:
    {
      /* C18O Photo-ionization, Photo-dissociation with shielding */
      theta = getcoss(NCO, NH2, 18);
      k = alpha * chi * exp (-gamm * av) * theta;
      break;
    }

    case 20:
      /* Depletion on the grains */
      {
        double prefreeze, n_gr;
        n_gr = grain_abundance * nh;
        prefreeze = M_PI * pow(grain_size,2.)*n_gr*sqrt(8.0*CONST_CGSM_BOLTZMANN/(M_PI*CONST_CGSM_MASS_PROTON));

        k = alpha*prefreeze*sqrt(tgas/beta);
        break;
      }

    case 21:
      /* Thermal desorption */
      {
        /* from Visser 2011 */
        double pregrain, n_gr, N_b;
        n_gr = grain_abundance * nh;
        N_b = GRAIN_SITES_PER_CM2 * M_PI * pow(grain_size,2.);
        // 2 layer desorption
        double f_spec = N_des/fmax(n_ice,N_des*N_b*n_gr);         
        pregrain = M_PI*pow(grain_size,2.)*n_gr*f_spec; 

        
        k = 4.0 * pregrain * alpha * exp(-gamm / tdust);
        

        break;
      }

    case 22:
      /* Cosmic ray desorption */
      {
        if (alpha == 0.0)
          {
            double v0 =
              pow (2 * GRAIN_SITES_PER_CM2 * gamm * CONST_CGSM_BOLTZMANN /
                   (M_PI * M_PI * beta * CONST_CGSM_MASS_PROTON),
                   0.5);
            k = v0 * FRACTION_TIME_GRAIN_70K * exp (-gamm / 70.);
            break;
          }
        else
          {
            k = alpha;
            break;
          }
      }

    case 23:
      /* Photo-desorption */
      {
	if (grain_abundance != 0)
	  {
      /* from Visser 2011 */
      double pregrain, n_gr, N_b;
      n_gr = grain_abundance * nh;
      N_b = GRAIN_SITES_PER_CM2 * M_PI * pow(grain_size,2.);
      // 2 layer desorption
      double f_spec = N_des/fmax(n_ice,N_des*N_b*n_gr);
      pregrain = M_PI*pow(grain_size,2.)*n_gr*f_spec;        

      // DRAINE_STANDARD_ISRF_FUV = 1.7e8 photons.cm^-2
      // COSMIC_DEFAULT = 1.3e-17 s^-1
      // 1e-4 to simulate CR desorption
      k = pregrain * alpha * DRAINE_STANDARD_ISRF_FUV * ((chi * exp(-3.02 * av)) + 1e-4*(cosx/5e-17));
	  }
	else
	  {
	    k = 0;
	  }
	break;
      }

    case 25:
      /* Hydrogenation on grains */
      {
        // fprintf(stdout, "n_hydro = %e\n", n_hydro);
        double prehydro, n_gr, N_b;
        n_gr = grain_abundance * nh;
        N_b = GRAIN_SITES_PER_CM2 * M_PI * pow(grain_size,2.);
        double f_prime = N_des/fmax(n_hydro,N_des*N_b*n_gr);
        // if (f_prime < MIN_FRAC) f_prime = 0.;
        // if (f_prime > 0.9999) f_prime = 1.;
        prehydro = M_PI*pow(grain_size,2.)*n_gr*f_prime;

        k = prehydro*sqrt(8.*CONST_CGSM_BOLTZMANN*tgas/(M_PI*CONST_CGSM_MASS_PROTON));
        break;
      }

    case 41:
      {
        k = alpha*(cosx/COSMIC_DEFAULT) * gamm/(1.0-0.5)*pow(tgas/300.0,beta);
        break;
      }

    case 42:
      /* cosmic ray CO photodissociation */
      {
        // calc k in solve.c (because we need to know amount of CO)
        k = 1;
        break;
      }

    case 43:
      {
        // calc k in solve.c (because we need to know amount of H,H2,He)
        k = 1;
        break;
      }

    case 60:
      {
        // Secondary ionization of H
        k = alpha*xray*0.56;
        break;
      }

    case 61:
      {
        // Secondary ionization of H2
        k = alpha*xray;
        break;
      }

    case 62:
      {
        //secondary ionization of other molecules
        k=alpha*xray;
        break;
      }

    case 70:
      {
        // photoelectron production from PAH
        // 1e-4 to simulate CR ionization
        k = alpha * (chi * exp(-3.02 * av) + 1.e-4);
        break;
      }

    case 71:
      {
        // Charge exchange with PAH
        k = alpha*0.5*pow(tgas/100.0,beta);
        break;
      }

    case 72:
      {
        // charge exchange with PAH
        k = alpha * 0.5 * pow(tgas/100.0,beta)*1.0/sqrt(gamm);
        break;
      }

    case 90:
    {
      // pumping of H2 --> H2*
      // get dissociation rate of H2 here, finish calc in solve.c
      double NH2_5e14, delta_v;
        NH2_5e14 = NH2 / 5.e14;
        delta_v = 2.e4; //thermal dispersion in cm/s
        theta = 0.965/pow(1.0 + (NH2_5e14/delta_v), 2.0) +
                  0.035/sqrt(1.0 + NH2_5e14) * exp(-8.5e-4*sqrt(1.0 + NH2_5e14));
        k = alpha * chi * exp (-gamm * av) * theta;
      break;
    }

    case 91:
    {
      // De-excitation H2* --> H2
      // get dissociation rate of H2 here, finish calc in solve.c
      double NH2_5e14, delta_v;
        NH2_5e14 = NH2 / 5.e14;
        delta_v = 2.e4; //thermal dispersion in cm/s
        theta = 0.965/pow(1.0 + (NH2_5e14/delta_v), 2.0) +
                  0.035/sqrt(1.0 + NH2_5e14) * exp(-8.5e-4*sqrt(1.0 + NH2_5e14));
        k = alpha * chi * exp (-gamm * av) * theta;
      break;
    }

    case 92:
    {
      // Reaction with H2*
      k = alpha * pow(tgas/300.0,beta)*exp(-fmax(0.0,gamm-30163.0)/tgas);
      break;
    }

    default:
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "unknown reaction type.\n");
      k = 0;
      break;
    }

  return k;
}







