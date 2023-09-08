/*
@author: Brian Kyanjo
@description: Solves a shallow water equation (swe) given a single left and right states
@date: 1st August 2023
@reference: Solver is described in J. Comput.Phys. (6): 3089-3113, March 2008 Augmented 
            Riemann Solvers for the swe with steady states and Inundation
*/

#ifndef GEOFLOOD_RIEMANN_UTILS_H
#define GEOFLOOD_RIEMANN_UTILS_H

#ifdef __cplusplus
extern "C" 
{
#if 0
}
#endif
#endif

extern __device__ void riemann_aug_JCP(int maxiter, int meqn, int mwaves, double hL, double hR,
                 double huL, double huR, double hvL, double hvR, double bL, double bR, 
                 double uL, double uR, double vL, double vR, double phiL, double phiR,
                 double sE1, double sE2, double drytol, double g, double* sw, double* fw);


extern __device__ void riemann_ssqfwave(int maxiter, int meqn, int mwaves, double hL, 
                 double hR, double huL, double huR, double hvL, double hvR, double bL, 
                 double bR, double uL, double uR, double vL, double vR, double phiL, 
                 double phiR, double sE1, double sE2, double drytol, double g, double* sw, double* fw) ;


extern __device__ void riemann_fwaves(int meqn, int mwaves, double hL, double hR, double huL, 
                double huR, double hvL, double hvR, double bL, double bR, double uL, double uR, 
                double vL, double vR, double phiL, double phiR, double s1, double s2, double drytol, double g, double* sw, double* fw) ;


extern __device__ void riemanntype(double hL, double hR, double uL, double uR, double hm, 
                 double s1m, double s2m, bool rare1, bool rare2, int maxiter,
                 double drytol, double g);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /*GEOFLOOD_RIEMANN_UTILS_H*/