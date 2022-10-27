/*
   solve.c - Build the ODE system and the jacobian matrix, and solve
   the system with CVODE.

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
#include <string.h>
#include <math.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#ifdef USE_LAPACK
#include <sunlinsol/sunlinsol_lapackdense.h>
#else
#include <sunlinsol/sunlinsol_dense.h>
#endif
#include <sundials/sundials_types.h>

#include "astrochem.h"
#include "network.h"
#include "rates.h"

static int f (realtype t, N_Vector y, N_Vector ydot, void *params);

static int jacobian (realtype t, N_Vector y, N_Vector fy,
                     SUNMatrix J, void *params, N_Vector tmp1,
                     N_Vector tmp2, N_Vector tmp3);

/*
  Right-hand side function of the ODE system.
*/

static int
f (realtype t __attribute__ ((unused)), N_Vector y, N_Vector ydot,
   void *params)
{
  int i;
  // fprintf(stdout, "%.5e\n", t);

  /* Get the parameters passed with *params. */

  double *reac_rates = ((params_t *) params)->reac_rates;
  const net_t *network = ((params_t *) params)->network;
  const react_t *reactions = ((params_t *) params)->reactions;
  int n_reactions = ((params_t *) params)->n_reactions;
  int n_species = ((params_t *) params)->n_species;
  double nh = ((params_t *) params)->nh;
  double av = ((params_t *) params)->av;
  double tgas = ((params_t *) params)->tgas;
  double tdust = ((params_t *) params)->tdust;
  double chi = ((params_t *) params)->chi;
  double cosmic = ((params_t *) params)->cosmic;
  double grain_size = ((params_t *) params)->grain_size;
  double grain_abundance = ((params_t *) params)->grain_abundance;
  double NCO = ((params_t *) params)->NCO;
  double NH2 = ((params_t *) params)->NH2;
  double NHD = ((params_t *) params)->NHD;
  double xray = ((params_t *) params)->xray;

  /* Loop over the reactions and build the right hand ODE system
     equations. */

  for (i = 0; i < n_species; i++)
    {
      NV_Ith_S (ydot, i) = 0;
    }

  /* get the hydrogenatable abundances */
  /* and ices */
  // FILE* fn;
  // FILE* ff;
  // FILE* fa;
  // FILE* fr;
  // FILE* fk;
  double n_hydro = 0.;
  double n_ice = 0.;
  for (i = 0; i < n_reactions; i++) {
    if (reactions[i].reaction_type == 25) {
      // if (NV_Ith_S (y, reactions[i].reactants[0]) > MIN_ICE_ABUND*nh) {
        n_hydro += NV_Ith_S (y, reactions[i].reactants[0]); // get the number of hydrogenatable ice
        // fa = fopen("abund.out","a");
        // fprintf(fa, "hydro %d\t%.5e\t%.5e\n", reactions[i].reaction_no,t,NV_Ith_S (y, reactions[i].reactants[0]));
        // fclose(fa);
      // }
    }
    if (reactions[i].reaction_type == 21) {
      // if (NV_Ith_S (y, reactions[i].reactants[0]) > MIN_ICE_ABUND*nh) {
        n_ice += NV_Ith_S (y, reactions[i].reactants[0]); // get the number of all ices
      // }
    }
    // fa = fopen("abund.out","a");
    // fprintf(fa, "ice %d\t%.5e\t%.5e\n", reactions[i].reaction_no,t,NV_Ith_S (y, reactions[i].reactants[0]));
    // fclose(fa);
  }
  
  // fn = fopen("n_ice.out","a");
  // fprintf(fn,"%.5e\t%.5e\t%.5e\t%.5e\n",t,n_ice,n_hydro,N_b*n_gr);
  // fclose(fn);


  for (i = 0; i < n_reactions; i++)
    {
      double y_product;

      /* Compute the production/destruction rate of the reaction. */

      y_product = 1;

      if (reactions[i].reaction_type == 0)
        {
          /* The formation of H2 is a first order reaction, so the
             product of the reactants needs to be divided by the H
             abundance. */
          y_product *= NV_Ith_S (y, reactions[i].reactants[0]);
          y_product *= reac_rates[i];
        }
      else if (reactions[i].reaction_type == 100)
        {
          /* HD formation */
          int ih = find_species("H",network);
          int id = find_species("D",network);
          double n_H, n_D;
          n_H = 0.;
          n_D = 0.;
          if (ih >= 0) n_H = NV_Ith_S(y,ih);
          if (id >= 0) n_D = NV_Ith_S(y,id);

          y_product *= reac_rates[i]/(1.e3 + fmax(n_H,n_D));
          y_product *= NV_Ith_S(y,reactions[i].reactants[0])*NV_Ith_S(y,reactions[i].reactants[1]);
        }
      else if (reactions[i].reaction_type == 21 || reactions[i].reaction_type == 23)
        {
          /* Desorption reactions depend on ice abundances */

          reac_rates[i] = rate (reactions[i].alpha,
                                reactions[i].beta,
                                reactions[i].gamma,
                                reactions[i].reaction_type,
                                reactions[i].reaction_no,
                                nh, NCO, NH2, NHD, av, tgas, tdust,
                                chi, cosmic,
                                grain_size,
                                grain_abundance,
                                NV_Ith_S (y, reactions[i].reactants[0]) / nh,
                                n_hydro, n_ice, xray);
          y_product *= reac_rates[i];
          double nJ = NV_Ith_S (y, reactions[i].reactants[0]);
          y_product *= nJ;
        }
      else if (reactions[i].reaction_type == 25)
      {
        /* hydrogenation also depends on ice thickness, so it is recomputed */
        reac_rates[i] = rate (reactions[i].alpha,
                                reactions[i].beta,
                                reactions[i].gamma,
                                reactions[i].reaction_type,
                                reactions[i].reaction_no,
                                nh, NCO, NH2, NHD, av, tgas, tdust,
                                chi, cosmic,
                                grain_size,
                                grain_abundance,
                                NV_Ith_S (y, reactions[i].reactants[0]) / nh,
                                n_hydro, n_ice, xray);
        y_product *= reac_rates[i];
        /* include n_H in rate */
        double nJ = NV_Ith_S (y, reactions[i].reactants[0]);
        y_product *= nJ;
        y_product *= NV_Ith_S (y, reactions[i].reactants[1]); //nH
      }
      else if (reactions[i].reaction_type == 42)
      {
        /* CR CO desorption depends on CO abundance */
        double co_abun, pd_eff;
        int ico = find_species("CO",network);
        if (ico < 0) { // no CO in network
          co_abun = 1.e-12;
        } else {
          co_abun = fabs(NV_Ith_S(y,ico)/nh)+1.e-12;
        }
        // sigmoid fit from Visser et al. (2018) to Heays et al. (2014)
        // n(CO)/nh + 1e-12 to avoid division by zero
        // always n(CO) to reflect CO self-shielding
        pd_eff = 56.14/(5.11e4*pow(co_abun,0.792)+1.)+4.3;
        reac_rates[i] = pd_eff*(cosmic+xray);

        y_product *= reac_rates[i];
        y_product *= NV_Ith_S (y, reactions[i].reactants[0]); // rate = k*n(CO{isotope})
      }
      else if (reactions[i].reaction_type == 43)
      {
        /* CR induced FUV reactions involving He processes */

        // from DALI code
          // from Staeuber/Doty code Eixion(3)/Eixion(1)
        int ihe = find_species("He",network);
        int ih2 = find_species("H2",network);
        int ih = find_species("H",network);
        double n_He = 0.;
        double n_H2 = 0.;
        double n_H = 0.;
        if (ihe >= 0) n_He = NV_Ith_S(y,ihe);
        if (ih2 >= 0) n_H2 = NV_Ith_S(y,ih2);
        if (ih >= 0) n_H = NV_Ith_S(y,ih);
        reac_rates[i] = reactions[i].alpha*((cosmic+xray)*0.0107)*n_He/
              (reactions[i].beta*n_H2+reactions[i].gamma*n_H+1e-20);

        y_product *= reac_rates[i];
        y_product *= NV_Ith_S (y, reactions[i].reactants[0]); // rate = k*n(j)

      }
      // else if (reactions[i].reaction_type == 90)
      // {
      //   /* pumping of H2 --> H2* */

      //   // get dissociation rate of H2
      //   // pass alpha, beta, gamm for H2 dissociation!
      //   reac_rates[i] = rate (5.18e-11,
      //                           1.00e+00,
      //                           3.02e+00,
      //                           reactions[i].reaction_type,
      //                           reactions[i].reaction_no,
      //                           nh, NCO, NH2, av, tgas, tdust, 
      //                           chi, cosmic,
      //                           grain_size,
      //                           grain_abundance,
      //                           NV_Ith_S (y, reactions[i].reactants[0]) / nh,
      //                           n_hydro, n_ice, xray);


      //   // see Visser et al. 2018, Le Bourlot et al. (1999)

      //   // get H,H2 number densities
      //   int ih = find_species("H",network);
      //   int ih2 = find_species("H2",network);
      //   double n_H = 0.;
      //   double n_H2 = 0.;
      //   if (ih2 >= 0) n_H2 = NV_Ith_S(y,ih2);
      //   if (ih >= 0) n_H = NV_Ith_S(y,ih);

      //   // downward collision rate
      //   double col_HH2,tinv;
      //   tinv = (tgas/1000.)+1.;
      //   col_HH2=pow(10.,-11.058+0.05554/tinv-2.3900/(tinv*tinv))*n_H + 
      //             pow(10.,-11.084-3.6706/tinv-2.0230/(tinv*tinv))*n_H2;
      //   reac_rates[i] *= 10.0;
      //   reac_rates[i] += col_HH2*exp(-30163.0/tgas);

      //   y_product *= reac_rates[i];
      //   y_product *= NV_Ith_S (y, reactions[i].reactants[0]);
      // }
      // else if (reactions[i].reaction_type == 91)
      // {
      //   /* de-excitation of H2* --> H2 */

      //   // get dissociation rate of H2
      //   // pass alpha, beta, gamm for H2 dissociation!
      //   reac_rates[i] = rate (5.18e-11,
      //                           1.00e+00,
      //                           3.02e+00,
      //                           reactions[i].reaction_type,
      //                           reactions[i].reaction_no,
      //                           nh, NCO, NH2, av, tgas, tdust, 
      //                           chi, cosmic,
      //                           grain_size,
      //                           grain_abundance,
      //                           NV_Ith_S (y, reactions[i].reactants[0]) / nh,
      //                           n_hydro, n_ice, xray);


      //   // see Visser et al. 2018, Le Bourlot et al. (1999)

      //   // get H,H2 number densities
      //   int ih = find_species("H",network);
      //   int ih2 = find_species("H2",network);
      //   double n_H = 0.;
      //   double n_H2 = 0.;
      //   if (ih2 >= 0) n_H2 = NV_Ith_S(y,ih2);
      //   if (ih >= 0) n_H = NV_Ith_S(y,ih);

      //   // downward collision rate
      //   double col_HH2,tinv;
      //   tinv = (tgas/1000.)+1.;
      //   col_HH2=pow(10.,-11.058+0.05554/tinv-2.3900/(tinv*tinv))*n_H + 
      //             pow(10.,-11.084-3.6706/tinv-2.0230/(tinv*tinv))*n_H2;
      //   reac_rates[i] *= 10.0;
      //   reac_rates[i] += 2.e-7+col_HH2;

      //   y_product *= reac_rates[i];
      //   y_product *= NV_Ith_S (y, reactions[i].reactants[0]);
      // }
      else
        {
          /* For other reactions, the production/destruction rate is
             the product of the reactants multiplied by the reaction rate. */
          int r;
          for( r = 0; r < MAX_REACTANTS; r++ )
            {
              if (reactions[i].reactants[r] != -1)
                y_product *= NV_Ith_S (y, reactions[i].reactants[r]);
            }
          y_product *= reac_rates[i];

        }
      /* Add a production term for each product of the reaction. */
      int p;
      for( p = 0 ; p < MAX_PRODUCTS; p++ )
        {
          if (reactions[i].products[p] != -1)
            NV_Ith_S (ydot, reactions[i].products[p]) += y_product;
        }
      /* Add a destruction term for each reactants of the reaction. */
      int r;
      for( r = 0; r < MAX_REACTANTS; r++ )
        {
          if (reactions[i].reactants[r] != -1)
            NV_Ith_S (ydot, reactions[i].reactants[r]) -= y_product;
        }
    }

  return (0);
}

/*
   Jacobian matrix.
   */

static int
jacobian (realtype t __attribute__ ((unused)), 
          N_Vector y,
          N_Vector fy __attribute__ ((unused)),
          SUNMatrix J, 
          void *params,
          N_Vector tmp1 __attribute__ ((unused)),
          N_Vector tmp2 __attribute__ ((unused)),
          N_Vector tmp3 __attribute__ ((unused)))
{
  int i, j;

  /* Get the parameters passed with *params. */

  double *reac_rates = ((params_t *) params)->reac_rates;
  const react_t *reactions = ((params_t *) params)->reactions;
  int n_reactions = ((params_t *) params)->n_reactions;
  int n_species = ((params_t *) params)->n_species;
  double nh = ((params_t *) params)->nh;
  double av = ((params_t *) params)->av;
  double chi = ((params_t *) params)->chi;
  double tdust = ((params_t *) params)->tdust;
  double tgas = ((params_t *) params)->tgas;
  double grain_size = ((params_t *) params)->grain_size;
  double grain_abundance = ((params_t *) params)->grain_abundance;
  double xray = ((params_t *) params)->xray;

  // FILE* fg;
  // fg = fopen("grain_abundance.out","a");
  // fprintf(fg, "%.5e\t%.5e\n", t,grain_abundance);
  // fclose(fg);

  double n_gr = grain_abundance * nh;
  double N_b = GRAIN_SITES_PER_CM2 * M_PI * pow(grain_size,2.);
  double N_des = 2.;

  // FILE* fj;

  double n_hydro = 0.;
  double n_ice = 0.;
  for (i = 0; i < n_reactions; i++) {
    if (reactions[i].reaction_type == 25) {
        n_hydro += NV_Ith_S (y, reactions[i].reactants[0]); // get the number of hydrogenatable ice
    }
    if (reactions[i].reaction_type == 21) {
        n_ice += NV_Ith_S (y, reactions[i].reactants[0]); // get the number of all ices
    }
  }

  /* Compute the jacobian matrix. */
  for (i = 0; i < n_species; i++)
    {
      for (j = 0; j < n_species; j++)
        {
          SM_ELEMENT_D (J, i, j) = 0;
        }
    }

  for (i = 0; i < n_reactions; i++)
    {

      /* Compute the Jacobian matrix elements corresponding to
         the reaction. */

      double y_product;

      if (reactions[i].reaction_type == 0)
        {
          /* H2 formation on grains */

          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -=
           2 * reac_rates[i];
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) +=
           reac_rates[i];
        }
      else if (reactions[i].reaction_type == 25) {
        /* hydrogenation */
        double dR_dn;
        double n_J = NV_Ith_S (y, reactions[i].reactants[0]);
        double n_H = NV_Ith_S (y, reactions[i].reactants[1]);
        if (n_hydro <= N_des*n_gr*N_b) {
          /* dR/dn_J   reactant 0 */
          dR_dn = reac_rates[i] * n_H;

          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].reactants[1], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) += dR_dn;

          /* dR/dn(H)  reactant 1 */
          dR_dn = reac_rates[i] * n_J;

          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[1]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].reactants[1], reactions[i].reactants[1]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[1]) += dR_dn;
        }
        else if (n_hydro > N_des*n_gr*N_b) {
          /* get f */
          double frac = n_J / n_hydro;

          /* dR/dn_J  reactant 0 */
          if (frac == 1.) {
            dR_dn = 0.;
          } else {
            dR_dn = reac_rates[i] * (1.-frac) * n_H;
          }

          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].reactants[1], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) += dR_dn;

          /* dR/dn(H)  reactant 1 */
          dR_dn = reac_rates[i] * n_J;

          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[1]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].reactants[1], reactions[i].reactants[1]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[1]) += dR_dn;
        }
      }

      else if (reactions[i].reaction_type == 23) {
        /* Photo-desorption */
        double dR_dn;
        double n_J = NV_Ith_S (y, reactions[i].reactants[0]);

        if (n_ice <= N_des*n_gr*N_b) {

          dR_dn = reac_rates[i];

          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) += dR_dn;
        }
        else if (n_ice > N_des*n_gr*N_b) {
          double frac = n_J / n_ice;

          if (frac == 1.) {
            dR_dn = 0.;
          } else {
            dR_dn = reac_rates[i] * (1.-frac);
          }
          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) += dR_dn;

        }
      }

      else if (reactions[i].reaction_type == 21) {
        /* thermal desorption */
        double dR_dn;
        double n_J = NV_Ith_S (y, reactions[i].reactants[0]);

        if (n_ice <= N_des*n_gr*N_b) {

          dR_dn = reac_rates[i];
          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) += dR_dn;


        }
        else if (n_ice > N_des*n_gr*N_b) {
          double frac = n_J / n_ice;

          if (frac == 1.) {
            dR_dn = 0.;
          } else {
            dR_dn = reac_rates[i] * (1.-frac);
          }
          SM_ELEMENT_D (J, reactions[i].reactants[0], reactions[i].reactants[0]) -= dR_dn;
          SM_ELEMENT_D (J, reactions[i].products[0], reactions[i].reactants[0]) += dR_dn;

        }
      }
      else
        {
          /* Other reactions */
          int r, r2, p;
          int nreactants = MAX_REACTANTS;
          for( r = 0; r < MAX_REACTANTS; r++ )
            {
              if( reactions[i].reactants[r] == -1 )
                {
                  nreactants = r;
                  break;
                }
            }
          for( r = 0; r < nreactants; r++ )
            {
              y_product =  reac_rates[i];
              for( r2 = 0; r2 < nreactants; r2++ )
                {
                  if ( r != r2 )
                    {
                      y_product *=  NV_Ith_S (y, reactions[i].reactants[ r2 ]);
                    }
                }

              for( r2 = 0; r2 < nreactants; r2++ )
                {
                  SM_ELEMENT_D (J, reactions[i].reactants[r2], reactions[i].reactants[r]) -= y_product;
                }
                // fj = fopen("jacob.out","a");
                // fprintf(fj, "%d\t%.5e\t%.5e\n", reactions[i].reaction_no,t,y_product);
                // fclose(fj);
              for( p = 0; p < MAX_PRODUCTS; p++ )
                {
                  if(  reactions[i].products[p] != -1 )
                    {
                      SM_ELEMENT_D (J, reactions[i].products[p], reactions[i].reactants[r]) += y_product;
                    }
                }
            }
        }
    }

  return (0);
}



/**
 * Initialize the solver
 */
int solver_init( const cell_t* cell, const net_t* network, const phys_t* phys,
                 const double* abundances , double density, double abs_err, double rel_err,
                 astrochem_mem_t* astrochem_mem )
{
  astrochem_mem->density = density;
  // fprintf(stdout, "cr ion = %.3e, xray ion = %.3e, total ion = %.3e\n", phys->cosmic, cell->xray, phys->cosmic+cell->xray);

  // FILE* fn;
  // FILE* ff;
  // FILE* fa;
  // FILE* fr;
  // FILE* fj;
  // FILE* fk;
  // fn = fopen("n_ice.out","w+");
  // ff = fopen("f_val.out","w+");
  // fa = fopen("abund.out","w+");
  // fr = fopen("rates.out","w+");
  // fj = fopen("jacob.out","w+");
  // fk = fopen("rateconst.out","w+");
  // fclose(fk);
  // fclose(fj);
  // fclose(fn);
  // fclose(ff);
  // fclose(fa);
  // fclose(fr);

  /* Allocate the work array and fill it with the initial
     concentrations. Ignore species that are not in the
     network. */

  if (( astrochem_mem->y = N_VNew_Serial (network->n_species)) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: %s\n", __FILE__, __LINE__,
               "array allocation failed.\n");
      return -1;
    }
  int i;
  for (i = 0; i < network->n_species; i++)
    {
      NV_Ith_S ( astrochem_mem->y, i) = abundances[i] * density; //this calculates the total number of elements for ODE solver
    }

  /* Allocate an array for the reaction rates and compute them. */

  if ((astrochem_mem->params.reac_rates =
       malloc (sizeof (double) * (unsigned int) network->n_reactions)) ==
      NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return -1;
    }


  for (i = 0; i < network->n_reactions; i++)
    {
      if (network->reactions[i].reaction_type == 90)
      {
        /* pumping of H2 --> H2* */

        // get dissociation rate of H2
        // pass alpha, beta, gamm for H2 dissociation!
        astrochem_mem->params.reac_rates[i] = rate (5.18e-11,
                                1.00e+00,
                                3.02e+00,
                                network->reactions[i].reaction_type,
                                network->reactions[i].reaction_no,
                                cell->nh, cell->NCO, cell->NH2, cell->NHD,
                                cell->av, cell->tgas, cell->tdust, 
                                phys->chi, phys->cosmic,
                                phys->grain_size,
                                phys->grain_abundance,
                                0., 0., 0., cell->xray);


        // see Visser et al. 2018, Le Bourlot et al. (1999)

        // get H,H2 number densities
        int ih = find_species("H",network);
        int ih2 = find_species("H2",network);
        double n_H = 0.;
        double n_H2 = 0.;
        if (ih2 >= 0) n_H2 = abundances[ih2]*cell->nh;
        if (ih >= 0) n_H = abundances[ih]*cell->nh;

        // downward collision rate
        double col_HH2,tinv;
        tinv = (cell->tgas/1000.)+1.;
        col_HH2=pow(10.,-11.058+0.05554/tinv-2.3900/(tinv*tinv))*n_H + 
                  pow(10.,-11.084-3.6706/tinv-2.0230/(tinv*tinv))*n_H2;
        astrochem_mem->params.reac_rates[i] *= 10.0;
        astrochem_mem->params.reac_rates[i] += col_HH2*exp(-30163.0/cell->tgas);
      }
      else if (network->reactions[i].reaction_type == 91)
      {
        /* de-excitation of H2* --> H2 */

        // get dissociation rate of H2
        // pass alpha, beta, gamm for H2 dissociation!
        astrochem_mem->params.reac_rates[i] = rate (5.18e-11,
                                1.00e+00,
                                3.02e+00,
                                network->reactions[i].reaction_type,
                                network->reactions[i].reaction_no,
                                cell->nh, cell->NCO, cell->NH2, cell->NHD,
                                cell->av, cell->tgas, cell->tdust, 
                                phys->chi, phys->cosmic,
                                phys->grain_size,
                                phys->grain_abundance,
                                0., 0., 0., cell->xray);


        // see Visser et al. 2018, Le Bourlot et al. (1999)

        // get H,H2 number densities
        int ih = find_species("H",network);
        int ih2 = find_species("H2",network);
        double n_H = 0.;
        double n_H2 = 0.;
        if (ih2 >= 0) n_H2 = abundances[ih2]*cell->nh;
        if (ih >= 0) n_H = abundances[ih]*cell->nh;

        // downward collision rate
        double col_HH2,tinv;
        tinv = (cell->tgas/1000.)+1.;
        col_HH2=pow(10.,-11.058+0.05554/tinv-2.3900/(tinv*tinv))*n_H + 
                  pow(10.,-11.084-3.6706/tinv-2.0230/(tinv*tinv))*n_H2;
        astrochem_mem->params.reac_rates[i] *= 10.0;
        astrochem_mem->params.reac_rates[i] += 2.e-7+col_HH2;
      }
      else {
      astrochem_mem->params.reac_rates[i] = rate (network->reactions[i].alpha,
                                                  network->reactions[i].beta,
                                                  network->reactions[i].gamma,
                                                  network->reactions[i].reaction_type,
                                                  network->reactions[i].reaction_no,
                                                  cell->nh, 
                                                  cell->NCO,
                                                  cell->NH2,
                                                  cell->NHD,
                                                  cell->av, 
                                                  cell->tgas,
                                                  cell->tdust, 
                                                  phys->chi,
                                                  phys->cosmic,
                                                  phys->grain_size,
                                                  phys->grain_abundance, 0., 0., 0., cell->xray);
      }
    }

  /* Fill out a structure containing the parameters of the function
     defining the ODE system and the jacobian. */

  astrochem_mem->params.network = network;
  astrochem_mem->params.reactions = network->reactions;
  astrochem_mem->params.n_reactions = network->n_reactions;
  astrochem_mem->params.n_species = network->n_species;
  astrochem_mem->params.nh = cell->nh;
  astrochem_mem->params.av = cell->av;
  astrochem_mem->params.tgas = cell->tgas;
  astrochem_mem->params.tdust = cell->tdust;
  astrochem_mem->params.NCO = cell->NCO;
  astrochem_mem->params.NH2 = cell->NH2;
  astrochem_mem->params.NHD = cell->NHD;
  astrochem_mem->params.chi = phys->chi;
  astrochem_mem->params.cosmic = phys->cosmic;
  astrochem_mem->params.grain_size = phys->grain_size;
  astrochem_mem->params.grain_abundance = phys->grain_abundance;
  astrochem_mem->params.xray = cell->xray;

  /* Define the ODE system and solve it using the Backward
     Differential Formulae method (BDF) with a Newton Iteration. The
     absolute error is multiplied by the density, because we compute
     concentrations and not abundances. */

  if (((astrochem_mem->cvode_mem = CVodeCreate (CV_BDF)) == NULL)
      || ((astrochem_mem->a = SUNDenseMatrix (network->n_species, network->n_species)) == NULL)
#ifdef USE_LAPACK
      || ((astrochem_mem->ls = SUNLinSol_LapackDense (astrochem_mem->y, astrochem_mem->a)) == NULL))
#else
      || ((astrochem_mem->ls = SUNLinSol_Dense (astrochem_mem->y, astrochem_mem->a)) == NULL))
#endif
    {
      fprintf (stderr, "astrochem: %s:%d: solver memory allocation failed.\n",
               __FILE__, __LINE__);
      return -1;
    }

  astrochem_mem->t = 0.0;
  if ((CVodeInit ( astrochem_mem->cvode_mem, f, astrochem_mem->t, astrochem_mem->y ) != CV_SUCCESS)
      || (CVodeSStolerances ( astrochem_mem->cvode_mem, rel_err, abs_err * density ) != CV_SUCCESS)
      || ((CVodeSetLinearSolver ( astrochem_mem->cvode_mem, astrochem_mem->ls, astrochem_mem->a) != CV_SUCCESS))
      || ((CVodeSetJacFn ( astrochem_mem->cvode_mem, jacobian) != CV_SUCCESS))
      || (CVodeSetUserData ( astrochem_mem->cvode_mem, &astrochem_mem->params) != CV_SUCCESS)
      || (CVodeSetMaxNumSteps (astrochem_mem->cvode_mem, CVODE_MXSTEPS) != CV_SUCCESS))
    {
      fprintf (stderr, "astrochem: %s:%d: solver initialization failed.\n",
               __FILE__, __LINE__);
      return -1;
    }
  return 0;
}

/*
  Solve the ODE system
 */

int solve( astrochem_mem_t* astrochem_mem, const net_t* network, double* abundances, double time , const cell_t* new_cell, int verbose )
{
  realtype t; /* Current time */

  /* Computing new reac rate and abundance if cell parameter have evolved
     since last call THIS IS ONLY RELEVENT FOR DYNAMIC*/
  if( new_cell != NULL )
    {
      int i,j;
      for (j = 0; j < network->n_species; j++)
        {
          NV_Ith_S (astrochem_mem->y, j) *= new_cell->nh / astrochem_mem->params.nh;
        }
        fprintf(stdout, "Setting new abundances in solve()\n");

      for (i = 0; i < network->n_reactions; i++)
        {
          astrochem_mem->params.reac_rates[i] = rate (network->reactions[i].alpha,
                                                      network->reactions[i].beta,
                                                      network->reactions[i].gamma,
                                                      network->reactions[i].reaction_type,
                                                      network->reactions[i].reaction_no,
                                                      new_cell->nh, new_cell->NCO, new_cell->NH2, new_cell->NHD,
                                                      new_cell->av, new_cell->tgas,
                                                      new_cell->tdust, astrochem_mem->params.chi,
                                                      astrochem_mem->params.cosmic,
                                                      astrochem_mem->params.grain_size,
                                                      astrochem_mem->params.grain_abundance, 0., 0., 0., new_cell->xray);
        }

      /* Fill out a structure containing the parameters of the function
         defining the ODE system and the jacobian. */

      astrochem_mem->params.nh = new_cell->nh;
      astrochem_mem->params.av = new_cell->av;
      astrochem_mem->params.tgas = new_cell->tgas;
      astrochem_mem->params.tdust = new_cell->tdust;
      astrochem_mem->params.NCO = new_cell->NCO;
      astrochem_mem->params.NH2 = new_cell->NH2;
      astrochem_mem->params.NHD = new_cell->NHD;
      astrochem_mem->params.xray = new_cell->xray;

      /* Re-initialize the solver */
      CVodeReInit (astrochem_mem->cvode_mem, astrochem_mem->t, astrochem_mem->y);
    }

  CVode ( astrochem_mem->cvode_mem, (realtype) time, astrochem_mem->y, &t, CV_NORMAL);
  astrochem_mem->t = time;
    
  /* Print the cell number, time and time step after each call. */

  if (verbose >= 2)
    {
      realtype h;

      CVodeGetLastStep (astrochem_mem->cvode_mem, &h);
      fprintf (stdout, "t = %8.2e  delta_t = %8.2e\n",
               (double) t / CONST_MKSA_YEAR,
               (double) h / CONST_MKSA_YEAR);
    }
  int i;
  for( i = 0; i < network->n_species; i++ )
    {
      abundances[i] =  (double) NV_Ith_S ( astrochem_mem->y , i ) / astrochem_mem->density ;
    }
  return 0;
}

/**
 * Close the solver
 */
void solver_close( astrochem_mem_t* astrochem_mem )
{
  if( astrochem_mem->params.reac_rates != NULL )
    {
      free (astrochem_mem->params.reac_rates);
    }
  if( astrochem_mem->y != NULL )
    {
      N_VDestroy_Serial (astrochem_mem->y);
    }
  if( astrochem_mem->cvode_mem != NULL )
    {
      CVodeFree (&astrochem_mem->cvode_mem);
    }
  if( astrochem_mem->ls != NULL )
    {
      SUNLinSolFree (astrochem_mem->ls);
    }
  if( astrochem_mem->a != NULL )
    {
      SUNMatDestroy (astrochem_mem->a);
    }
}

/**
 *  Allocate the abundance vector
 */
int alloc_abundances ( const net_t* network, double** abundances )
{
  if (( *abundances =
        malloc (sizeof (double) * network->n_species )) == NULL)
    {
      fprintf (stderr, "astrochem: %s:%d: array allocation failed.\n",
               __FILE__, __LINE__);
      return -1;
    }
  int i;
  for( i = 0 ; i < network->n_species ; i++ )
    {
      (*abundances)[i]=0;
    }
  return 0;
}

/**
 *  Initialise specific abundances
 */
int set_initial_abundances( const char** species, int n_initialized_species, const double* initial_abundances,
                            const net_t* network, double* abundances )
{
  int i,j;
  for( i = 0 ; i < network->n_alloc_species ; i++ )
    {
      for( j = 0 ; j < n_initialized_species ; j++ )
        {
          if( strcmp( network->species[i].name , species[j] ) == 0 )
            {
              abundances[i] = initial_abundances[j];
            }
        }
    }
  return 0;
}

/**
 * Free the abundances
 */
void free_abundances ( double* abundances )
{
  if( abundances != NULL )
    {
      free ( abundances);
    }
}
