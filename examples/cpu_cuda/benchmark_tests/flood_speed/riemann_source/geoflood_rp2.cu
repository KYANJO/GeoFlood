/*
@author: Brian Kyanjo
@description: Approximate Riemann Solver for the Shallow Water Equations in Cuda
@date: November 28th 2023
@reference: Solver is described in J. Comput.Phys. (6): 3089-3113, March 2008 Augmented 
            Riemann Solvers for the swe with steady states and Inundation
*/

#define maxiter 1

#include "../flood_speed_user.h"
#include <math.h>
#include <fc2d_cudaclaw.h>
#include <fc2d_cudaclaw_check.h>
#include <fc2d_cudaclaw_options.h>
#include <cudaclaw_user_fort.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_include_all.h>

/* === Begin fuction riemann_aug_JCP======================================================== @description: - Solves swe give single left and right states
@note: - To use the original solver call with maxiter=1.
       - This solver allows iteration when maxiter > 1. The iteration seems to help  
         with instabilities that arise (with any solver) as flow becomes transcritical 
         over variable topography due to loss of hyperbolicity. 
*/

__device__ void riemann_aug_JCP(int maxiter, int meqn, int mwaves, double hL,
    double hR, double huL, double huR, double hvL, double hvR, 
    double bL, double bR, double uL, double uR, double vL, 
    double vR, double phiL, double phiR, double sE1, double sE2, 
    double drytol, double g, double* sw, double* fw)
{


}

/* === Begin fuction Riemann type ============
 @description: Determines the Riemann structure (wave-type in each family)
*/

__device__ void riemanntype(double hL, double hR, double uL, double uR, double *hm, 
                            double *s1m, double *s2m, bool *rare1, bool *rare2, int maxiter,
                            double drytol, double g)
{
    // Local variables
    double u1m, u2m, um, h0, F_max, F_min, dfdh, F0, slope, gL, gR;
    double sqrtgh1, sqrtgh2;
    double h0;
    int iter; 

    // Test for Riemann structure
    double h_min = fmin(hR,hL);
    double h_max = fmax(hR,hL);
    double delu = uR - uL;

    /* Have dry state on either side 
    - Only one rarefaction wave
    - another shock wave has 0 jump and moves at the same speed as one edge of the    rarefaction  wave */
    if (h_min <= drytol)
    {
        *hm = 0.0;
        um = 0.0;
    
        /* Either hR or hL is almost zero, so the expression below corresponds
           to either Eqn. (54a) or Eqn. (54b) in the JCP paper */
        *s1m = uR + uL - 2.0 * sqrt(g * hR) + 2.0 * sqrt(g * hL);
        *s2m = *s1m; 
        *rare1 = (hL <= 0.0) ? false : true;
        *rare2 = !(*rare1);
    } else {
        F_min = delu + 2.0 * (sqrt(g * h_min) - sqrt(g * h_max));
        F_max = delu + (h_max - h_min) * sqrt(0.5 * g * (h_max + h_min) / (h_max * h_min));

        if (F_min > 0.0){  // 2-rarefactions
            /* Eqn (13.56) in the FVMHP book */
            *hm = (1.0 / (16.0 * g)) * pow(fmax(0.0, -delu + 2.0 * (sqrt(g * hL) + sqrt(g * hR))), 2);
            um = copysign(1.0, *hm) * (uL + 2.0 * (sqrt(g * hL) - sqrt(g * *hm)));
            *s1m = uL + 2.0 * sqrt(g * hL) - 3.0 * sqrt(g * *hm);
            *s2m = uR - 2.0 * sqrt(g * hR) + 3.0 * sqrt(g * *hm);
            *rare1 = true;
            *rare2 = true;
        } else if (F_max <= 0.0) { // 2-shocks
            /* Below it solves for the intersection of two Hugoniot loci to get the
            accurate Riemann solution */
            /* Root finding using a Newton iteration on sqrt(h) */
            h0 = h_max;
            for (iter = 1; iter <= MAXITER; iter++) {
                gL = sqrt(0.5 * g * (1.0 / h0 + 1.0 / hL));
                gR = sqrt(0.5 * g * (1.0 / h0 + 1.0 / hR));
                F0 = delu + (h0 - hL) * gL + (h0 - hR) * gR;
                dfdh = gL - g * (h0 - hL) / (4.0 * h0 * h0 * gL) + gR - g * (h0 - hR) / (4.0 * h0 * h0 * gR);
                slope = 2.0 * sqrt(h0) * dfdh;
                h0 = pow(sqrt(h0) - F0 / slope, 2);
            }
            *hm = h0;
            /* u1m and u2m are Eqns (13.19) and (13.20) in the FVMHP book */
            u1m = uL - (*hm - hL) * sqrt(0.5 * g * (1.0 / *hm + 1.0 / hL));
            u2m = uR + (*hm - hR) * sqrt(0.5 * g * (1.0 / *hm + 1.0 / hR));
            um = 0.5 * (u1m + u2m);
            *s1m = u1m - sqrt(g * *hm);
            *s2m = u2m + sqrt(g * *hm);
            *rare1 = false;
            *rare2 = false;
        } else { // 1-shock or 1-rarefaction
            h0 = h_min;
            for (iter = 1; iter <= MAXITER; iter++) {
                F0 = delu + 2.0 * (sqrt(g * h0) - sqrt(g * h_max)) + (h0 - h_min) * sqrt(0.5 * g * (1.0 / h0 + 1.0 / h_min));
                slope = (F_max - F0) / (h_max - h_min);
                h0 = h0 - F0 / slope;
            }
            *hm = h0;
            sqrtgh2 = sqrt(g * *hm);
            if (hL > hR) {
                sqrtgh1 = sqrt(g * hL);
                /* Eqn (13.55) in the FVMHP book */
                um = uL + 2.0 * sqrtgh1 - 2.0 * sqrtgh2;
                *s1m = uL + 2.0 * sqrtgh1 - 3.0 * sqrtgh2;
                *s2m = uL + 2.0 * sqrtgh1 - sqrtgh2;

                *rare1 = true;
                *rare2 = false;
            } else {
                sqrtgh1 = sqrt(g * hR);
                um = uR - 2.0 * sqrtgh1 + 2.0 * sqrtgh2;
                *s1m = uR - 2.0 * sqrtgh1 + sqrtgh2;
                *s2m = uR - 2.0 * sqrtgh1 + 3.0 * sqrtgh2;
                *rare1 = false;
                *rare2 = true;
            }
        }
    }
} /* End of riemanntype function */
