/* 
@author: David L. George
@modified to C by: Brian Kyanjo
@date: 31 July 2023
@description: Solves normal Riemann problems for the 2D shallow water equations (swe) with 
topography:
            h_t + (hu)_x + (hv)_y = 0
            (hu)_t + (hu^2 + 1/2gh^2)_x + (huv)_y = -ghb_x
            (hv)_t + (huv)_x + (hv^2 + 1/2gh^2)_y = -ghb_y
where h is the height, u is the x velocity, v is the y velocity, g is the gravitational constant, and b is the topography.
@input: ql - conatins the state vector at the left edge of each cell
        qr - contains the state vector at the right edge of each cell
        
        This data is along a slice in the x-direction if idir = 0 or along a slice in the y-direction if idir = 1.

        idir - indicates the direction of the slice

@note: - The ith Riemann problem has left state qr(i-1,:) and right state ql(i,:).
       - This solver allows the user to easily select a Riemann solver in riemann_solvers.c,    this routine initializes all the variables for the swe, accounting for wet dry boundary, dry cells, wave speeds, etc.
       
@reference: David L. George
*/

#define maxiter 1

#include "../flood_speed_user.h"
#include "variables.h"
#include <math.h>
#include <fc2d_geoclaw.h>
#include <fc2d_cudaclaw_check.h>
#include <fc2d_cudaclaw_options.h>
#include <cudaclaw_user_fort.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_include_all.h>

/* Extern declarations*/
extern __constant__ GeofloodVars d_geofloodVars;

/* function prototypes */
__device__ void riemanntype(double hL, double hR, double uL, double uR, double *hm, double *s1m, double *s2m, bool *rare1, bool *rare2);

__device__ void riemann_type(double hL, double hR, double uL, double uR, double hm, 
    double s1m, double s2m, bool rare1, bool rare2);

__device__ void riemann_aug_JCP(int meqn, int mwaves, double hL,
    double hR, double huL, double huR, double hvL, double hvR, 
    double bL, double bR, double uL, double uR, double vL, 
    double vR, double phiL, double phiR, double sE1, double sE2, double* sw, double* fw, int ix, int iy, int idir);

__device__ void riemann_aug_JCP(int meqn, int mwaves, double hL,
        double hR, double huL, double huR, double hvL, double hvR, 
        double bL, double bR, double uL, double uR, double vL, 
        double vR, double phiL, double phiR, double sE1, double sE2, double* sw1, double* sw2, double* sw3, double* fw11, double* fw12, double* fw13, double* fw21, double* fw22, double* fw23, double* fw31, double* fw32, double* fw33, int ix, int iy, int idir);

__device__ void flood_speed_compute_speeds(int idir, int meqn, int mwaves, int maux,
                                            double ql[], double  qr[],
                                            double auxl[], double auxr[],
                                            double s[])
{

    /* Access the __constant__ variables in variables.h */
    double s_grav = d_geofloodVars.gravity;
    double drytol = d_geofloodVars.dry_tolerance;
    double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = d_geofloodVars.coordinate_system;
    int mcapa = d_geofloodVars.mcapa;

    int mu = 1+idir;
    // int mv = 2-idir;

    double hhat = (ql[0] + qr[0])/2.0;
    double hsq2 = sqrt(ql[0]) + sqrt(qr[0]);
    double uhat = (ql[mu]/sqrt(ql[0]) + qr[mu]/sqrt(qr[0]))/hsq2;
    // double vhat = (ql[mv]/sqrt(ql[0]) + qr[mv]/sqrt(qr[0]))/hsq2;
    double chat = sqrt(s_grav*hhat);

    // Roe wave speeds
    double roe1 = uhat - chat;
    double roe3 = uhat + chat;

    // left and right state wave speeds
    double s1l = ql[mu]/ql[0] - sqrt(s_grav*ql[0]);
    double s3r = qr[mu]/qr[0] + sqrt(s_grav*qr[0]);

    // Einfeldt wave speeds
    double s1 = fmin(s1l,roe1);
    double s3 = fmax(s3r,roe3);

    double s2 = 0.5*(s1+s3);

    s[0] = s1;
    s[1] = s2;
    s[2] = s3;
}

__device__ cudaclaw_cuda_speeds_t flood_speed_speeds = flood_speed_compute_speeds;

void flood_speed_assign_speeds(cudaclaw_cuda_speeds_t *speeds)
{
    cudaError_t ce = cudaMemcpyFromSymbol(speeds, flood_speed_speeds, sizeof(cudaclaw_cuda_speeds_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (flood_speed_compute_speeds): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}

/* Normal Riemann solver for the 2d shallow water equations with topography */
// __device__ void cudaflood_rpn2(int idir, int meqn, int mwaves,
//                                 int maux, double ql[], double qr[],
//                                 double auxl[], double auxr[],
//                                 double fwave[], double s[], 
//                                 double amdq[], double apdq[], int ix, int iy)
__device__ void cudaflood_rpn2(int idir, int meqn, int mwaves,
                                int maux, double ql[], double qr[],
                                double auxl[], double auxr[],
                                double fwave[], double s[], 
                                double amd_q[], double apd_q[], int ix, int iy)
{
    /* Access the __constant__ variables in variables.h */
    double s_grav = d_geofloodVars.gravity;
    double drytol = d_geofloodVars.dry_tolerance;
    double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = d_geofloodVars.coordinate_system;
    int mcapa = d_geofloodVars.mcapa; 

    /* Local variables */
    double wall[3], fw[9], sw[3];
    double hR, hL, huR, huL, hvR, hvL, uR, uL, vR, vL, phiR, phiL;
    double bR, bL, sL, sR, sRoe1, sRoe2, sE1, sE2, uhat, chat;
    double hstar, hstartest, dxdc;
    double s1m, s2m;
    bool rare1, rare2;
    int mw, mu, mv;

    /* vectorized version */
    double fw11, fw12, fw13, fw21, fw22, fw23, fw31, fw32, fw33;
    double sw1, sw2, sw3;

    /* Swapping left to right  (cudaclaw_flux2.cu)*/
    // double *qr = q_l;
    // double *ql = q_r;
    // double *auxr = aux_l;
    // double *auxl = aux_r;
    double *amdq = apd_q;
    double *apdq = amd_q;

       //   print at only one thread
    //    int mx = 16, my = 16, mbc = 2;
    // //    int thread_index = threadIdx.x;
    //    int ifaces_x, ifaces_y;
    //   ifaces_x = mx + 2*mbc-1;
    //   ifaces_y = my + 2*mbc-1;
    //   int ix = thread_index % ifaces_x;
    //   int iy = thread_index/ifaces_x;  

      bool debug;
      if (idir == 0)
      {
        debug = 1;
      }
      else{
        debug = 0;
      }

    /* === Initializing === */
    /* inform of a bad riemann problem from the start */
    if ((qr[0] < 0.0) || (ql[0] < 0.0)) {
        printf("Negative input: hl, hr = %f,%f\n", ql[0], qr[0]);
    }

    // Initialize Riemann problem for the grid interface 
    for (mw=0; mw<mwaves; ++mw)
    {
        s[mw] = 0.0;
        fwave[mw + 0*mwaves] = 0.0;
        fwave[mw + 1*mwaves] = 0.0; 
        fwave[mw + 2*mwaves] = 0.0;
    }

    /* set normal direction */
    mu = 1+idir;
    mv = 2-idir;

    /* zero (small) negative values if they exist */
    // left state
    if (qr[0] < 0.0) {
        qr[0] = 0.0;
        qr[1] = 0.0;
        qr[2] = 0.0;
    }

    // right state
    if (ql[0] < 0.0) {
        ql[0] = 0.0;
        ql[1] = 0.0;
        ql[2] = 0.0;
    }

    // if (debug){
    //     printf("ix = %d, iy = %d\n " \ 
    //     "qr[0] = %.16f, ql[0] = %.16f\n" \
    //     "qr[1] = %.16f, ql[1] = %.16f\n" \
    //     "qr[2] = %.16f, ql[2] = %.16f\n\n", ix,iy,qr[0],ql[0],qr[1],ql[1],qr[2],ql[2]);
    // }

    // // Skip problem if in a completely dry area
    if (qr[0] <= drytol && ql[0] <= drytol) {
        goto label30;
    }

    // if (ql[0] > drytol || qr[0] > drytol) {
        /* Riemann problem variables */
        // hL  = qr[0];
        // hR  = ql[0];
        // huL = qr[mu];
        // huR = ql[mu];
        // bL = auxr[0];
        // bR = auxl[0];
        hL = ql[0];
        hR = qr[0];
        huL = ql[mu];
        huR = qr[mu];
        bL = auxl[0];
        bR = auxr[0];

        // hvL = qr[mv];
        // hvR = ql[mv];
        hvL = ql[mv];
        hvR = qr[mv];

        // Check for wet/dry left boundary
        if (hR > drytol) {
            uR = huR / hR;
            vR = hvR / hR;
            phiR = 0.5 * s_grav * (hR * hR) + (huR * huR) / hR;
        } else {
            hR = 0.0;
            huR = 0.0;
            hvR = 0.0;
            uR = 0.0;
            vR = 0.0;
            phiR = 0.0;
        }

        // Check for wet/dry right boundary
        if (hL > drytol) {
            uL = huL / hL;
            vL = hvL / hL;
            phiL = 0.5 * s_grav * (hL * hL) + (huL * huL) / hL;
        } else {
            hL = 0.0;
            huL = 0.0;
            hvL = 0.0;
            uL = 0.0;
            vL = 0.0;
            phiL = 0.0;
        }
    
        // if (debug){
        //     printf("ix = %d, iy = %d\n " \ 
        //     "hL = %.16f, hR = %.16f\n" \
        //     "huL = %.16f, huR = %.16f\n" \
        //     "hvL = %.16f, hvR = %.16f\n" \
        //     "uL = %.16f, uR = %.16f\n" \
        //     "vL = %.16f, vR = %.16f\n" \
        //     "phiL = %.16f, phiR = %.16f\n" \
        //     "bL = %.16f, bR = %.16f\n\n", ix,iy,hL,hR,huL,huR,hvL,hvR,uL,uR,vL,vR,phiL,phiR,bL,bR);
        // }

        /* left and right surfaces depth inrelation to topography */
        wall[0] = 1.0;
        wall[1] = 1.0;
        wall[2] = 1.0;
        if (hR <= drytol) {
            /* determine the wave structure */
            riemanntype(hL, hL, uL, -uL, &hstar, &s1m, &s2m, &rare1, &rare2);
            // riemann_type(hL, hL, -uL, uL, hstar, s1m, s2m, rare1, rare2);

        //     if (debug){
        //     printf("ix = %d, iy = %d\n " \ 
        //     "hL = %.16f, uL = %.16f\n" \
        //     "hstar = %.16f\n" \
        //     "s1m = %.16f, s2m = %.16f\n" \
        //     "rare1 = %d, rare2 = %d\n\n", ix,iy,hL,uL,hstar,s1m,s2m,rare1,rare2);
        // }

            hstartest = fmax(hL,hstar);
            if (hstartest + bL < bR) {
                /* hL+bL < bR and hstar+bL < bR, so water can't overtop right cell 
                (move into right cell) so right state should become ghost values 
                that mirror left for wall problem) */
                wall[1] = 0.0;
                wall[2] = 0.0;
                hR = hL;
                huR = -huL;
                bR = bL;
                phiR = phiL;
                uR = -uL;
                vR = vL;
                /* here we already have huR =- huL, so we don't need to change it */
            } else if (hL+bL < bR) {
                /* hL+bL < bR and hstar+bL >bR, so we set bR to the water level in 
                the left cell so that water can possibly overtop the right cell (move into the right cell) */ 
                bR = hL + bL;
            }
        } else if (hL <= drytol) { /* right surface is lower than left topo */
            /* determine the Riemann structure */
            riemanntype(hR, hR, -uR, uR, &hstar, &s1m, &s2m, &rare1, &rare2);
            // riemann_type(hR, hR, uR, -uR, hstar, s1m, s2m, rare1, rare2);
            hstartest = fmax(hR,hstar);

            // if (debug){
            //     printf("ix = %d, iy = %d\n " \ 
            //     "hR = %.16f, uR = %.16f\n" \
            //     "hstar = %.16f\n" \
            //     "s1m = %.16f, s2m = %.16f\n" \
            //     "rare1 = %d, rare2 = %d\n\n", ix,iy,hR,uR,hstar,s1m,s2m,rare1,rare2);
            // }

            if (hstartest + bR < bL) //left state should become ghost values that mirror right for wall problem
            {
                wall[0] = 0.0;
                wall[1] = 0.0;
                hL = hR;
                huL = -huR;
                bL = bR;
                phiL = phiR;
                uL = -uR;
                vL = vR;
            } else if (hR+bR < bL) {
                bL = hR + bR;
            }
        }

        // if (debug){
        //     printf("ix = %d, iy = %d\n " \ 
        //     "hL = %.16f, hR = %.16f\n" \
        //     "huL = %.16f, huR = %.16f\n" \
        //     "hvL = %.16f, hvR = %.16f\n" \
        //     "uL = %.16f, uR = %.16f\n" \
        //     "vL = %.16f, vR = %.16f\n" \
        //     "phiL = %.16f, phiR = %.16f\n" \
        //     "bL = %.16f, bR = %.16f\n\n", ix,iy,hL,hR,huL,huR,hvL,hvR,uL,uR,vL,vR,phiL,phiR,bL,bR);
        // }

        /* determine wave speeds */
        sL = uL - sqrt(s_grav*hL); // 1 wave speed of left state
        sR = uR + sqrt(s_grav*hR); // 2 wave speed of right state

        uhat = (sqrt(s_grav*hL)*uL + sqrt(s_grav*hR)*uR)/(sqrt(s_grav*hL) + sqrt(s_grav*hR)); // Roe average
        chat = sqrt(0.5*s_grav*(hL+hR)); // Roe average
        sRoe1 = uhat - chat; // Roe wave speed 1 wave
        sRoe2 = uhat + chat; // Roe wave speed 2 wave

        sE1 = fmin(sL,sRoe1); // Einfeldt wave speed 1 wave
        sE2 = fmax(sR,sRoe2); // Einfeldt wave speed 2 wave

        // if (ix == 7 && iy == 15) {
        //     if (debug){
        //         printf("ix = %d, iy = %d\n " \ 
        //         "sL = %.16f, sR = %.16f\n" \ 
        //         "sRoe1 = %.16f, sRoe2 = %.16f\n" \ 
        //         "sE1 = %.16f, sE2 = %.16f\n\n", ix,iy,sL,sR,sRoe1,sRoe2,sE1,sE2);
        //     }
        // }
        /* --- end of initializing --- */

        /* === solve Riemann problem === */
        // riemann_aug_JCP(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,sw,fw,ix,iy,idir);
        riemann_aug_JCP(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,&sw1,&sw2,&sw3,&fw11,&fw12,&fw13,&fw21,&fw22,&fw23,&fw31,&fw32,&fw33,ix,iy,idir);
        
        // Debugging check for NaNs 
        // if (tid == 0) {
        //     printf("hL = %e, hR = %e\n", hL, hR);
        //     printf("huL = %e, huR = %e\n", huL, huR);
        //     printf("hvL = %e, hvR = %e\n", hvL, hvR);
        //     printf("uL = %e, uR = %e\n", uL, uR);
        //     printf("vL = %e, vR = %e\n", vL, vR);
        //     printf("phiL = %e, phiR = %e\n", phiL, phiR);
        //     printf("bL = %e, bR = %e\n", bL, bR);
        // }


        // eliminate ghost fluxes for wall
        // for (mw=0; mw<3; mw++) {
        //     sw[mw] *= wall[mw];
        //     fw[mw] *= wall[mw];
        //     fw[mw + 1*mwaves] *= wall[mw];
        //     fw[mw + 2*mwaves] *= wall[mw];
        // }
        sw1 *= wall[0];
        sw2 *= wall[1];
        sw3 *= wall[2];
        fw11 *= wall[0];
        fw21 *= wall[0];
        fw31 *= wall[0];
        fw12 *= wall[1];
        fw22 *= wall[1];
        fw32 *= wall[1];
        fw13 *= wall[2];
        fw23 *= wall[2];
        fw33 *= wall[2];    
        
        // update fwave and corresponding speeds
        // for (mw=0; mw<mwaves; mw++) {
        //     s[mw] = sw[mw];
        //     fwave[mw] = fw[mw];
        //     fwave[mw + mu*mwaves] = fw[mw + 1*mwaves];
        //     fwave[mw + mv*mwaves] = fw[mw + 2*mwaves];
        // }
        s[0] = sw1;
        s[1] = sw2;
        s[2] = sw3;
        fwave[0] = fw11;
        fwave[1] = fw21;
        fwave[2] = fw31;
        fwave[3] = fw12;
        fwave[4] = fw22;
        fwave[5] = fw32;
        fwave[6] = fw13;
        fwave[7] = fw23;
        fwave[8] = fw33;

    // }

    // Debugging check for NaNs
    // if (tid == 0) {
    //     printf("s1 = %e, s2 = %e, s3 = %e\n", s[0], s[1], s[2]);
    //     printf("fwave[0] = %e, fwave[1] = %e, fwave[2] = %e\n", fwave[0], fwave[1], fwave[2]);
    //     printf("fwave[3] = %e, fwave[4] = %e, fwave[5] = %e\n", fwave[3], fwave[4], fwave[5]);
    //     printf("fwave[6] = %e, fwave[7] = %e, fwave[8] = %e\n", fwave[6], fwave[7], fwave[8]);
    // }

    label30: // (similar to 30 continue in Fortran)

    /* --- Capacity or Mapping from Latitude Longitude to physical space ----*/
    if (mcapa > 0) {
        if (idir == 0) {
            dxdc = earth_radius*deg2rad;
        } else {
            dxdc = earth_radius*cos(auxl[2])*deg2rad;
        }

        // update fwave and corresponding speeds
        for (mw=0; mw<mwaves; mw++) {
            s[mw] = dxdc*s[mw];
            fwave[mw] = dxdc*fwave[mw];
            fwave[mw + mwaves] = dxdc*fwave[mw + mwaves];
            fwave[mw + 2*mwaves] = dxdc*fwave[mw + 2*mwaves];
        }
    }

    /* --- compute fluctuations --- */
    for (mw=0; mw<mwaves; mw++) {
        if (s[mw] < 0.0) {
            amdq[mw] += fwave[mw];
            amdq[mw] += fwave[mw + 1*mwaves];
            amdq[mw] += fwave[mw + 2*mwaves];
        } else if (s[mw] > 0.0) {
            apdq[mw] += fwave[mw];
            apdq[mw] += fwave[mw + 1*mwaves];
            apdq[mw] += fwave[mw + 2*mwaves];
        } else {
            amdq[mw] += 0.5*fwave[mw];
            amdq[mw] += 0.5*fwave[mw + 1*mwaves];
            amdq[mw] += 0.5*fwave[mw + 2*mwaves];
            apdq[mw] += 0.5*fwave[mw];
            apdq[mw] += 0.5*fwave[mw + 1*mwaves];
            apdq[mw] += 0.5*fwave[mw + 2*mwaves];
        }
    }
    // Debugging check for NaNs
    // if (tid == 0) {
    //     printf("amdq[0] = %e, amdq[1] = %e, amdq[2] = %e\n", amdq[0], amdq[1], amdq[2]);
    //     printf("apdq[0] = %e, apdq[1] = %e, apdq[2] = %e\n", apdq[0], apdq[1], apdq[2]);
    // }
}


__device__ cudaclaw_cuda_rpn2_t flood_speed_rpn2 = cudaflood_rpn2;

void flood_speed_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, flood_speed_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (cudaflood_rpn2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}

/* Transverse Riemann solver for the 2d shallow water equations with topography 
@desc: Using The Jacobian matrix from left cell (imp == 0) or right cell (imp == 1) to compute the transverse fluxes.
*/

__device__ void cudaflood_rpt2(int idir, int meqn, int mwaves, int maux,
                double q_l[], double q_r[], double aux1[], 
                double aux2[], double aux3[], int imp, 
                double asdq[], double bmasd_q[], double bpasd_q[], int ix, int iy) 
{
    /* Access the __constant__ variables in variables.h */
    double s_grav = d_geofloodVars.gravity;
    double drytol = d_geofloodVars.dry_tolerance;
    double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = d_geofloodVars.coordinate_system;
    int mcapa = d_geofloodVars.mcapa;

    int mw, mu, mv;
    double s[3], beta[3];
    double r[3][3];
    double h, u, v;
    double delf1, delf2, delf3;
    double dxdcm, dxdcp, topo1, topo3, eta;

    /* Swapping left to right  (cudaclaw_flux2.cu)*/
    double *qr = q_l;
    double *ql = q_r;
    double *bmasdq = bpasd_q;
    double *bpasdq = bmasd_q;

    mu = 1+idir;
    mv = 2-idir;

    h = (imp == 0) ? qr[0] : ql[0];

    bool debug = (idir == 0) ? 1 : 0;
  
    // if (h <= drytol) return; // skip problem if dry cell (leaves bmadsq(:) = bpasdq(:) = 0)
    if (h > drytol) {
        /* Compute velocities in relevant cell, and other quantities */
        if (imp == 0) {
            // fluctuations being split is left-going
            u = qr[mu] / h;
            v = qr[mv] / h;
            eta = h + aux2[0];
            topo1 = aux1[0];
            topo3 = aux3[0];
        } else {
            // fluctuations being split is right-going
            u = ql[mu] / h;
            v = ql[mv] / h;
            eta = h + aux2[0];
            topo1 = aux1[0];
            topo3 = aux3[0];
        }

        // Debugging check for NaNs
        // if (tid == 0) {
        //     printf("h = %e, u = %e, v = %e, eta = %e, topo1 = %e, topo3 = %e\n", h, u, v, eta, topo1, topo3);
        // }

        /* Check if cell that transverse wave go into are both too high: */
        // if (eta < fmin(topo1, topo3)) return; 
        if (eta >= fmin(topo1, topo3)) {

            /* Check if cell that transverse waves go into are both to high, if so,
            do the splitting (no dry cells), and compute necessary quantities */
            if (coordinate_system == 2) {
                // On the sphere
                if (idir == 1) {
                    dxdcp = earth_radius * deg2rad;
                    dxdcm = dxdcp;
                } else {
                    if (imp == 0) {
                        dxdcp = earth_radius * cos(aux3[2]) * deg2rad;
                        dxdcm = earth_radius * cos(aux1[2]) * deg2rad;
                    } else {
                        dxdcp = earth_radius * cos(aux3[2]) * deg2rad;
                        dxdcm = earth_radius * cos(aux1[2]) * deg2rad;
                    }
                }
            } else {
                // Cartesian
                dxdcp = 1.0;
                dxdcm = 1.0;
            }

            /* Compute some speeds necessary for the Jacobian 
            - Computing upgoing, downgoing waves either in cell on left (if imp==0)
                or on the right (if imp==1) 
            - To achieve this we use q values in cells above and below, however these
                aren't available (only in aux values)
            */
            s[0] = v - sqrt(s_grav * h);
            s[1] = v;
            s[2] = v + sqrt(s_grav * h);

            // Debugging check for NaNs
            if (ix == 7 && iy == 15) {
                // if ((hL >= 0.3280909317849093) && (hR >= 0.3280909317849093)){
                if (debug){
                    printf("ix = %d, iy = %d\n " \
                    "s[0] = %e, s[1] = %e, s[2] = %e\n", ix,iy, s[0], s[1], s[2]);
                }
            }

            /* Determine asdq decomposition (beta) */
            delf1 = asdq[0];
            delf2 = asdq[mu];
            delf3 = asdq[mv];

            // Debugging check for NaNs
            // if (debug) {
            //     printf("delf1 = %e, delf2 = %e, delf3 = %e\n", delf1, delf2, delf3);
            // }

            beta[0] = (s[2]*delf1 - delf3) / (s[2] - s[0]);
            beta[1] = -u*delf1 + delf2;
            beta[2] = (delf3 - s[0]*delf1) / (s[2] - s[0]);

            /* set-up eigenvectors */
            // r[0] = 1.0;
            // r[mu] = u;
            // r[mv] = s[0];
            r[0][0] = 1.0;
            r[1][0] = u;
            r[2][0] = s[0];

            // r[0 + mwaves] = 0.0;
            // r[mu + mwaves] = 1.0;
            // r[mv + mwaves] = 0.0;
            r[0][1] = 0.0;
            r[1][1] = 1.0;
            r[2][1] = 0.0;

            // r[0 + 2*mwaves] = 1.0;
            // r[mu + 2*mwaves] = u;
            // r[mv + 2*mwaves] = s[2];
            r[0][2] = 1.0;
            r[1][2] = u;
            r[2][2] = s[2];

            // Debugging check for NaNs
            // if (tid == 0) {

            //     printf("beta[0] = %e, beta[1] = %e, beta[2] = %e\n", beta[0], beta[1], beta[2]);
            // }

            /* Compute transverse fluctuations */
            for (mw = 0; mw < 3; mw++) {
                if ((s[mw] < 0.0) && (eta >= topo1)) {
                    // bmasdq[0] += dxdcm * s[mw]*beta[mw]*r[mw + mwaves];
                    // bmasdq[mu] += dxdcm * s[mw]*beta[mw]*r[mw + mwaves];
                    // bmasdq[mv] += dxdcm * s[mw]*beta[mw]*r[mw + 2*mwaves];
                    bmasdq[0] += dxdcm * s[mw]*beta[mw]*r[0][mw];
                    bmasdq[mu] += dxdcm * s[mw]*beta[mw]*r[1][mw];
                    bmasdq[mv] += dxdcm * s[mw]*beta[mw]*r[2][mw];
                } else if ((s[mw] > 0.0) && (eta >= topo3)) {
                    // bpasdq[0] += dxdcp * s[mw]*beta[mw]*r[mw + mwaves];
                    // bpasdq[mu] += dxdcp * s[mw]*beta[mw]*r[mw + mwaves];
                    // bpasdq[mv] += dxdcp * s[mw]*beta[mw]*r[mw + 2*mwaves];
                    bpasdq[0] += dxdcp * s[mw]*beta[mw]*r[0][mw];
                    bpasdq[mu] += dxdcp * s[mw]*beta[mw]*r[1][mw];
                    bpasdq[mv] += dxdcp * s[mw]*beta[mw]*r[2][mw];
                }
            }
            // Debugging check for NaNs
            // if (tid == 0) {
            //     printf("bmasdq[0] = %e, bmasdq[1] = %e, bmasdq[2] = %e\n", bmasdq[0], bmasdq[1], bmasdq[2]);
            //     printf("bpasdq[0] = %e, bpasdq[1] = %e, bpasdq[2] = %e\n", bpasdq[0], bpasdq[1], bpasdq[2]);
            // }
        }
    }
}



__device__ cudaclaw_cuda_rpt2_t flood_speed_rpt2 = cudaflood_rpt2;

void flood_speed_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, flood_speed_rpt2, sizeof(cudaclaw_cuda_rpt2_t));

    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (cudaflood_rpt2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}

/* === Begin fuction riemann_aug_JCP======================================================== @description: - Solves swe give single left and right states
@note: - To use the original solver call with maxiter=1.
       - This solver allows iteration when maxiter > 1. The iteration seems to help  
         with instabilities that arise (with any solver) as flow becomes transcritical 
         over variable topography due to loss of hyperbolicity. 
*/

// __device__ void riemann_aug_JCP(int meqn, int mwaves, double hL,
//     double hR, double huL, double huR, double hvL, double hvR, 
//     double bL, double bR, double uL, double uR, double vL, 
//     double vR, double phiL, double phiR, double sE1, double sE2, double* sw, double* fw, int ix, int iy, int idir)
__device__ void riemann_aug_JCP(int meqn, int mwaves, double hL,
        double hR, double huL, double huR, double hvL, double hvR, 
        double bL, double bR, double uL, double uR, double vL, 
        double vR, double phiL, double phiR, double sE1, double sE2, double* sw1, double* sw2, double* sw3, double* fw11, double* fw12, double* fw13, double* fw21, double* fw22, double* fw23, double* fw31, double* fw32, double* fw33, int ix, int iy, int idir)
{
    /* Access the __constant__ variables in variables.h */
    double s_grav = d_geofloodVars.gravity;
    double drytol = d_geofloodVars.dry_tolerance;
    double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = d_geofloodVars.coordinate_system;
    int mcapa = d_geofloodVars.mcapa;

    /* Local variables */
    // double A[9], r[9], lambda[3], del[3], beta[3];
    double lambda[3], beta[3],del[3];
    double A[3][3], r[3][3];
    double delh, delhu, delphi, delb, delnorm;
    double rare1st, rare2st, sdelta, raremin, raremax;
    double criticaltol, convergencetol;
    double criticaltol_2, hustar_interface;
    double s1s2bar, s1s2tilde, hbar, hLstar, hRstar;
    double huRstar, huLstar, uRstar, uLstar, hstarHLL;
    double deldelh, deldelphi;
    double s1m, s2m, hm;
    double det1, det2, det3, determinant;
    bool rare1, rare2, rarecorrector, rarecorrectortest, sonic;
    int mw, k, iter;

    int mu = 1; // x-direction
    int mv = 2; // y-direction

    bool debug = (idir == 0) ? 1 : 0;

    /* determine del vectors */
    delh = hR - hL;
    delhu = huR - huL;
    delphi = phiR - phiL;
    delb = bR - bL;
    delnorm = delh * delh + delphi * delphi;

    /* Determine the Riemann structure */
    riemanntype(hL,hR,uL,uR,&hm,&s1m,&s2m,&rare1,&rare2);
    // riemann_type(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2);

    if (ix == 7 && iy == 15) {
        // if ((hL >= 0.3280909317849093) && (hR >= 0.3280909317849093)){
            if (debug){
                printf("ix = %d, iy = %d\n " \ 
                "hL = %.16f, hR = %.16f\n" \
                "uL = %.16f, uR = %.16f\n" \
                "hm = %.16f\n" \
                "s1m = %.16f, s2m = %.16f\n" \
                "g = %.16f, drytol = %.16f\n" \
                "rare1 = %d, rare2 = %d\n\n", ix,iy,hL,hR,uL,uR,hm,s1m,s2m,s_grav,drytol,rare1,rare2);
            }
        // }
    }
   
    /* For the solver to handle depth negativity, depth dh is included in the decompostion which gives as acess to using the depth positive semidefinite solver (HLLE). This makes the system to have 3 waves instead of 2. where the 1st and 3rd are the eigenpairs are related to the flux Jacobian matrix of the original SWE (since s1<s2<s3, and have been modified by Einfeldt to handle depth non-negativity) and the 2nd is refered to as the the entropy corrector wave since its introduced to correct entropy violating solutions with only 2 waves. */
    
    /* The 1st and 3rd speeds are the eigenvalues of the Jacobian matrix of the original SWE modified by Einfeldt's for use with the HLLE solver. */
    lambda[0] = fmin(sE1, s2m); /* sE1 - flux Jacobian eigen value s2m - Roe speed */
    lambda[2] = fmax(sE2, s1m); /* sE2 - flux Jacobian eigen value s1m - Roe speed*/

    /* Einfeldt's speeds */
    sE1 = lambda[0]; 
    sE2 = lambda[2];

    /* The 2nd speed is the entropy corrector wave speed. */
    lambda[1] = 0.0; /* no strong or significant rarefaction waves */
    
    /* determine the middle state in the HLLE solver */
    hstarHLL = fmax((hL*uL - hR*uR + (sE2 * hR) - (sE1 * hL)) / (sE2 - sE1), 0.0); /* middle state between the two discontinuities (positive semidefinite depth) */

    /* === determine the middle entropy corrector wave === */
    /* rarecorrectortest = .true. provides a more accurate Riemann solution but is more expensive. This is because a nonlinear Riemann solution with  2 nonlinear waves as a linear Riemann solution 3 (or 2 jump discontionuities to approximate 1 smooth nonlinear rarefaction if it's large). When rarecorrectortest = .false. the approximate solution has only 2 jump discontinuities instead of 3, so its less accurate but faster. */
    rarecorrectortest = false;
    rarecorrector = false;
    if (rarecorrectortest) {
        sdelta = lambda[2] - lambda[0];
        raremin = 0.5; /* indicate a large rarefaction wave but not large */
        raremax = 0.9; /* indicate a very large rarefaction wave */
       /* i.e (the total speed difference between the fastest and slowest wave in the Riemann solution = 0.5) */

        if (rare1 && sE1 * s1m < 0.0) raremin = 0.2;
        if (rare2 && sE2 * s2m < 0.0) raremin = 0.2;

        if (rare1 || rare2) {
            /* check which rarefaction is the strongest */
            rare1st = 3.0 * (sqrt(s_grav * hL) - sqrt(s_grav * hm));
            rare2st = 3.0 * (sqrt(s_grav * hR) - sqrt(s_grav * hm));
            if (fmax(rare1st, rare2st) > raremin * sdelta && fmax(rare1st, rare2st) < raremax * sdelta) {
                rarecorrector = true;
                if (rare1st > rare2st) {
                    lambda[1] = s1m;
                } else if (rare2st > rare1st) {
                    lambda[1] = s2m;
                } else {
                    lambda[1] = 0.5 * (s1m + s2m);
                }
            }
        }
        if (hstarHLL < fmin(hL, hR) / 5.0) rarecorrector = false;
    }

    /* determining modified eigen vectors */
    for (mw = 0; mw < mwaves; mw++) {   
        // r[mw] = 1.0; 
        // r[mw + mwaves] = lambda[mw]; 
        // r[mw + 2*mwaves] = pow(lambda[mw],2.0);
        r[0][mw] = 1.0;
        r[1][mw] = lambda[mw];
        r[2][mw] = pow(lambda[mw],2.0);
    }

    /* no strong rarefaction wave */
    if (!rarecorrector) {
        lambda[1]= 0.5*(lambda[0] + lambda[2]);
        // r[mwaves] = 0.0; // r[0,1]
        // r[mwaves + mu] = 0.0; // r[1,1]
        // r[mwaves + mv] = 1.0; // r[2,1]
        r[0][1] = 0.0;
        r[1][1] = 0.0;
        r[2][1] = 1.0;
    }

    /* === Determine the steady state wave === */
    criticaltol = fmax(drytol*s_grav, 1.0e-6);
    criticaltol_2 = sqrt(criticaltol);
    deldelh = -delb;
    deldelphi = -0.5 * (hR + hL) * (s_grav * delb); /* some approximation of the source term \int_{x_{l}}^{x_{r}} -g h b_x dx */

    /* determine a few quantities needed for steady state wave if iterated */
    hLstar = hL;
    hRstar = hR;
    uLstar = uL;
    uRstar = uR;
    huLstar = uLstar * hLstar;
    huRstar = uRstar * hRstar;

    /* iterate to better find the steady state wave */
    convergencetol = 1e-6;
    for (iter=1; iter <= maxiter; iter++) {
        /* determine steady state wave (this will be subtracted from the delta vectors */
        if (fmin(hLstar,hRstar) < drytol && rarecorrector) {
            rarecorrector = false;
            hLstar = hL;
            hRstar = hR;
            uLstar = uL;
            uRstar = uR;
            huLstar = uLstar*hLstar;
            huRstar = uRstar*hRstar;
            lambda[1] = 0.5*(lambda[0] + lambda[2]);
            // r[mwaves] = 0.0; // r[0,1]
            // r[mwaves + mu] = 0.0; // r[1,1]
            // r[mwaves + mv] = 1.0; // r[2,1]
            r[0][1] = 0.0;
            r[1][1] = 0.0;
            r[2][1] = 1.0;
        }

        /* For any two states; Q_i and Q_i-1, eigen values of SWE must satify: lambda(q_i)*lambda(q_i-1) = u^2 -gh, writing this conditon as a function of Q_i and Q_i-1, u and h become averages in lambda(q_i)*lambda(q_i-1) = u^2 -gh and these averages are denoted by bar and tilde. */
        hbar = fmax(0.5 * (hLstar + hRstar), 0.0);
        s1s2bar = 0.25 * pow((uLstar + uRstar),2) - (s_grav * hbar);
        s1s2tilde = fmax(0.0, uLstar * uRstar) - (s_grav * hbar);

        /* Based on the above conditon, smooth staedy state over slopping bathymetry cannot have a sonic point. Therefore, for regions with monotonically varying bathymetry, steady-state flow is either entirely subsonic (-u^2 +gh > 0) or entirely supersonic. */
        sonic = false;
        if (fabs(s1s2bar) <= criticaltol) {
            sonic = true;
        } else if (s1s2bar * s1s2tilde <= criticaltol * criticaltol) {
            sonic = true;
        } else if (s1s2bar * sE1 * sE2 <= criticaltol * criticaltol) {
            sonic = true;
        } else if (fmin(fabs(sE1), fabs(sE2)) < criticaltol_2) {
            sonic = true;
        } else if (sE1 < criticaltol_2 && s1m > -criticaltol_2) {
            sonic = true;
        } else if (sE2 > -criticaltol_2 && s2m < criticaltol_2) {
            sonic = true;
        } else if ((uL + sqrt(s_grav * hL)) * (uR + sqrt(s_grav * hR)) < 0.0) {
            sonic = true;
        } else if ((uL - sqrt(s_grav * hL)) * (uR - sqrt(s_grav * hR)) < 0.0) {
            sonic = true;
        }

        /* find jump in h, deldelh */
        if (sonic) {
            deldelh = -delb;
        } else {
            deldelh = delb * s_grav * hbar / s1s2bar;
        }

        /* find bounds in case of critical state resonance, or negative states */
        if (sE1 < -criticaltol && sE2 > criticaltol) {
            deldelh = fmin(deldelh, hstarHLL * (sE2 - sE1) / sE2);
            deldelh = fmax(deldelh, hstarHLL * (sE2 - sE1) / sE1);
        } else if (sE1 >= criticaltol) {
            deldelh = fmin(deldelh, hstarHLL * (sE2 - sE1) / sE1);
            deldelh = fmax(deldelh, -hL);
        } else if (sE2 <= -criticaltol) {
            deldelh = fmin(deldelh, hR);
            deldelh = fmax(deldelh, hstarHLL * (sE2 - sE1) / sE2);
        }

        /* find jump in phi, ddphi */
        if (sonic) {
            deldelphi = -s_grav * hbar * delb;
        } else {
            deldelphi = -delb * s_grav * hbar * s1s2tilde / s1s2bar;
        }

        /* find bounds in case of critical state resonance, or negative states */
        deldelphi = fmin(deldelphi, s_grav * fmax(-hLstar * delb, -hRstar * delb));
        deldelphi = fmax(deldelphi, s_grav * fmin(-hLstar * delb, -hRstar * delb));

        /* determine the delta vectors */
        del[0] = delh - deldelh;
        del[1] = delhu;
        del[2] = delphi - deldelphi;  

        /* Determine coefficients beta(k) using crammer's rule
          first determine the determinant of the eigenvector matrix */
        // det1 = r[0]*(r[mwaves + mu]*r[2*mwaves + mv] - r[2*mwaves + mu]*r[mwaves + mv]);
        // det2 = r[mwaves]*(r[mu]*r[2*mwaves + mv] - r[2*mwaves + mu]*r[mv]);
        // det3 = r[2*mwaves]*(r[mu]*r[mwaves + mv] - r[mwaves + mu]*r[mv]);
        det1 = r[0][0]*(r[1][1]*r[2][2] - r[1][2]*r[2][1]);
        det2 = r[0][1]*(r[1][0]*r[2][2] - r[1][2]*r[2][0]);
        det3 = r[0][2]*(r[1][0]*r[2][1] - r[1][1]*r[2][0]);
        determinant = det1 - det2 + det3;

        /* solve for beta(k) */
        for(k=0; k < 3; k++)
        {   
            for(mw=0; mw < 3; mw++)
            {
                // A[mw] = r[mw]; 
                // A[mw + mwaves] = r[mw + mwaves];
                // A[mw + 2*mwaves] = r[mw + 2*mwaves];
                A[0][mw] = r[0][mw];
                A[1][mw] = r[1][mw];
                A[2][mw] = r[2][mw];
            }
            // A[k] = del[0];
            // A[mwaves + k] = del[1];
            // A[2*mwaves + k] = del[2];
            // det1 = A[0]*(A[mwaves + mu]*A[2*mwaves + mv] - A[2*mwaves + mu]*A[mwaves + mv]);
            // det2 = A[mwaves]*(A[mu]*A[2*mwaves + mv] - A[2*mwaves + mu]*A[mv]);
            // det3 = A[2*mwaves]*(A[mu]*A[mwaves + mv] - A[mwaves + mu]*A[mv]);
            A[0][k] = del[0];
            A[1][k] = del[1];
            A[2][k] = del[2];
            det1 = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]);
            det2 = A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
            det3 = A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
            beta[k] = (det1 - det2 + det3)/determinant;
        }

        /* exit if things aren't changing */
        if (fabs(pow(del[0],2)+pow(del[2],2.0) - delnorm) < convergencetol) break;

        delnorm = pow(del[0],2)+pow(del[2],2.0); /* update delnorm */

        /* find new states qLstar and qRstar on either side of interface */
        hLstar = hL;
        hRstar = hR;
        uLstar = uL;
        uRstar = uR;
        huLstar = uLstar*hLstar;
        huRstar = uRstar*hRstar;

        /* left state depth and momentum updates */
        for (mw=0; mw < mwaves; mw++)
        {
            if (lambda[mw] < 0.0)
            {
                // hLstar = hLstar + beta[mw]*r[mw]; 
                // huLstar = huLstar + beta[mw]*r[mw + mwaves]; 
               hLstar = hLstar + beta[mw]*r[0][mw];
               huLstar = huLstar + beta[mw]*r[1][mw];
            }
        }

        /* right state depth and momentum updates */
        for (mw = mwaves-1; mw >= 0; mw--)
        {
            if (lambda[mw] > 0.0)
            {
                // hRstar = hRstar - beta[mw]*r[mw]; 
                // huRstar = huRstar - beta[mw]*r[mw + mwaves]; 
                hRstar = hRstar - beta[mw]*r[0][mw];
                huRstar = huRstar - beta[mw]*r[1][mw];
            }
        }

        /* left state velocity update */
        if (hLstar > drytol) 
        {
            uLstar = huLstar/hLstar;
        }
        else  /* dry state */
        {
            hLstar = fmax(hLstar,0.0);
            uLstar = 0.0;
        }

        /* right state velocity update */
        if (hRstar > drytol) 
        {
            uRstar = huRstar/hRstar;
        }
        else /* dry state */
        {
            hRstar = fmax(hRstar,0.0);
            uRstar = 0.0;
        }
    } /* end of  iteration on the Riemann problem*/

    /* === determine the fwaves and speeds=== */
    // for (mw=0; mw < mwaves; mw++)
    // {
    //     sw[mw] = lambda[mw];
    //     // fw[mw] = beta[mw]*r[mw + mwaves]; 
    //     // fw[mw + mwaves] = beta[mw]*r[mw + 2*mwaves]; 
    //     // fw[mw + 2*mwaves] = beta[mw]*r[mw + mwaves]; 
    //     fw[mw] = beta[mw]*r[1][mw];
    //     fw[mw + mwaves] = beta[mw]*r[2][mw];
    //     fw[mw + 2*mwaves] = beta[mw]*r[1][mw];
    // }
    *sw1 = lambda[0];
    *sw2 = lambda[1];
    *sw3 = lambda[2];
    // mw = 0;
    *fw11 = beta[0]*r[1][0];
    *fw21 = beta[0]*r[2][0];
    *fw31 = beta[0]*r[1][0];
    // mw = 1;
    *fw12 = beta[1]*r[1][1];
    *fw22 = beta[1]*r[2][1];
    *fw32 = beta[1]*r[1][1];
    // mw = 2;
    *fw13 = beta[2]*r[1][2];
    *fw23 = beta[2]*r[2][2];
    *fw33 = beta[2]*r[1][2];

    // find transverse components (ie huv jumps)
    // fw[mv] *= vL;
    // fw[2*mwaves + mv] *= vR;
    // fw[mwaves + mv] = 0.0;
    *fw31 *= vL;
    *fw33 *= vR;
    *fw32 = 0.0;

    // hustar_interface = hL*uL + fw[0];
    // if (hustar_interface <= 0.0) {
    //     fw[mv] += (hR * uR * vR - hL * uL * vL - fw[mv] - fw[2*mwaves + mv]);
    // } else {
    //     fw[2*mwaves + mv] += (hR * uR * vR - hL * uL * vL - fw[mv] - fw[2*mwaves + mv]);
    // }
    hustar_interface = hL*uL + *fw11;
    if (hustar_interface <= 0.0) {
        *fw31 += (hR * uR * vR - hL * uL * vL - *fw31 - *fw33);
    } else {
        *fw33 += (hR * uR * vR - hL * uL * vL - *fw31 - *fw33);
    }

}

/* === Begin fuction Riemann type ============
 @description: Determines the Riemann structure (wave-type in each family)
*/

__device__ void riemanntype(double hL, double hR, double uL, double uR, double *hm, 
                            double *s1m, double *s2m, bool *rare1, bool *rare2)
{
    /* Access the __constant__ variables in variables.h */
    double s_grav = d_geofloodVars.gravity;
    double drytol = d_geofloodVars.dry_tolerance;
    double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = d_geofloodVars.coordinate_system;
    int mcapa = d_geofloodVars.mcapa;

    // Local variables
    double um, u1m, u2m, h0, F_max, F_min, dfdh, F0, slope, gL, gR;
    double sqrtgh1, sqrtgh2;
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
        *s1m = uR + uL - 2.0 * sqrt(s_grav * hR) + 2.0 * sqrt(s_grav * hL);
        *s2m = *s1m; 
        *rare1 = (hL <= 0.0) ? false : true;
        *rare2 = !(*rare1);
    } else {
        F_min = delu + 2.0 * (sqrt(s_grav * h_min) - sqrt(s_grav * h_max));
        F_max = delu + (h_max - h_min) * sqrt(0.5 * s_grav * (h_max + h_min) / (h_max * h_min));

        if (F_min > 0.0){  // 2-rarefactions
            /* Eqn (13.56) in the FVMHP book */
            *hm = (1.0 / (16.0 * s_grav)) * pow(fmax(0.0, -delu + 2.0 * (sqrt(s_grav * hL) + sqrt(s_grav * hR))), 2);
            um = copysign(1.0, *hm) * (uL + 2.0 * (sqrt(s_grav * hL) - sqrt(s_grav * *hm)));
            *s1m = uL + 2.0 * sqrt(s_grav * hL) - 3.0 * sqrt(s_grav * *hm);
            *s2m = uR - 2.0 * sqrt(s_grav * hR) + 3.0 * sqrt(s_grav * *hm);
            *rare1 = true;
            *rare2 = true;
        } else if (F_max <= 0.0) { // 2-shocks
            /* Below it solves for the intersection of two Hugoniot loci to get the
            accurate Riemann solution */
            /* Root finding using a Newton iteration on sqrt(h) */
            h0 = h_max;
            for (iter = 1; iter <= maxiter; iter++) {
                gL = sqrt(0.5 * s_grav * (1.0 / h0 + 1.0 / hL));
                gR = sqrt(0.5 * s_grav * (1.0 / h0 + 1.0 / hR));
                F0 = delu + (h0 - hL) * gL + (h0 - hR) * gR;
                dfdh = gL - s_grav * (h0 - hL) / (4.0 * h0 * h0 * gL) + gR - s_grav * (h0 - hR) / (4.0 * h0 * h0 * gR);
                slope = 2.0 * sqrt(h0) * dfdh;
                h0 = pow(sqrt(h0) - F0 / slope, 2);
            }
            *hm = h0;
            /* u1m and u2m are Eqns (13.19) and (13.20) in the FVMHP book */
            u1m = uL - (*hm - hL) * sqrt(0.5 * s_grav * (1.0 / *hm + 1.0 / hL));
            u2m = uR + (*hm - hR) * sqrt(0.5 * s_grav * (1.0 / *hm + 1.0 / hR));
            um = 0.5 * (u1m + u2m);
            *s1m = u1m - sqrt(s_grav * *hm);
            *s2m = u2m + sqrt(s_grav * *hm);
            *rare1 = false;
            *rare2 = false;
        } else { // 1-shock or 1-rarefaction
            h0 = h_min;
            for (iter = 1; iter <= maxiter; iter++) {
                F0 = delu + 2.0 * (sqrt(s_grav * h0) - sqrt(s_grav * h_max)) + (h0 - h_min) * sqrt(0.5 * s_grav * (1.0 / h0 + 1.0 / h_min));
                slope = (F_max - F0) / (h_max - h_min);
                h0 = h0 - F0 / slope;
            }
            *hm = h0;
            sqrtgh2 = sqrt(s_grav * *hm);
            if (hL > hR) {
                sqrtgh1 = sqrt(s_grav * hL);
                /* Eqn (13.55) in the FVMHP book */
                um = uL + 2.0 * sqrtgh1 - 2.0 * sqrtgh2;
                *s1m = uL + 2.0 * sqrtgh1 - 3.0 * sqrtgh2;
                *s2m = uL + 2.0 * sqrtgh1 - sqrtgh2;

                *rare1 = true;
                *rare2 = false;
            } else {
                sqrtgh1 = sqrt(s_grav * hR);
                um = uR - 2.0 * sqrtgh1 + 2.0 * sqrtgh2;
                *s1m = uR - 2.0 * sqrtgh1 + sqrtgh2;
                *s2m = uR - 2.0 * sqrtgh1 + 3.0 * sqrtgh2;
                *rare1 = false;
                *rare2 = true;
            }
        }
    }
} /* End of riemanntype function */



__device__  void riemann_type(double hL, double hR, double uL, double uR, double hm, 
    double s1m, double s2m, bool rare1, bool rare2)
{
    double g = d_geofloodVars.gravity;
    double drytol = d_geofloodVars.dry_tolerance;

    double um,u1m,u2m,delu;
    double h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR;
    double sqrtgh1,sqrtgh2;
    int iter;

    /* Test for Riemann structure */
    h_min = min(hR,hL);
    h_max = max(hR,hL);
    delu = uR - uL;

   if (h_min <= drytol){
    hm = 0.0;
    um = 0.0;
    s1m = uR + uL - 2.0*sqrt(g*hR) + 2.0*sqrt(g*hL);
    s2m = uR + uL - 2.0*sqrt(g*hR) + 2.0*sqrt(g*hL);

    if (hL <= 0.0){
        rare1 = false;
        rare2 = true;
    } else {
        rare1 = true;
        rare2 = false;
    }
   } else {
        F_min = delu + 2.0*(sqrt(g*h_min) - sqrt(g*h_max));
        F_max = delu + (h_max - h_min)*(sqrt(0.5*g*(h_max + h_min)/(h_max*h_min)));

        if (F_min < 0.0) {
            /* 2-rarefactions */
            hm = (1.0/(16.0*g))*(pow(max(0.0,-delu + 2.0*(sqrt(g*hL) + sqrt(g*hR))),2));
            um = copysign(1.0,hm)*(uL + 2.0*(sqrt(g*hL) - sqrt(g*hm)));
            s1m = uL + 2.0*sqrt(g*hL) - 3.0*sqrt(g*hm);
            s2m = uR - 2.0*sqrt(g*hR) + 3.0*sqrt(g*hm);
            rare1 = true;
            rare2 = true;
        } else if (F_max <= 0.0) {
            /* 2-shocks */
            /* Root finding using a Newton iteration on sqrt(h) */
            h0 = h_max;
            for (iter = 1; iter <= maxiter; iter++) {
                gL = sqrt(0.5*g*(1.0/h0 + 1.0/hL));
                gR = sqrt(0.5*g*(1.0/h0 + 1.0/hR));
                F0 = delu + (h0 - hL)*gL + (h0 - hR)*gR;
                dfdh = gL - g*(h0 - hL)/(4.0*h0*h0*gL) + gR - g*(h0 - hR)/(4.0*h0*h0*gR);
                slope = 2.0*sqrt(h0)*dfdh;
                h0 = pow(sqrt(h0) - F0/slope,2);
            }
            hm = h0;
            /* u1m and u2m are Eqns (13.19) and (13.20) in the FVMHP book */
            u1m = uL - (hm - hL)*sqrt(0.5*g*(1.0/hm + 1.0/hL));
            u2m = uR + (hm - hR)*sqrt(0.5*g*(1.0/hm + 1.0/hR));
            um = 0.5*(u1m + u2m);
            s1m = u1m - sqrt(g*hm);
            s2m = u2m + sqrt(g*hm);
            rare1 = false;
            rare2 = false;
        } else {
            /* 1-shock or 1-rarefaction */
            h0 = h_min;
            for (iter = 1; iter <= maxiter; iter++) {
                F0 = delu + 2.0*(sqrt(g*h0) - sqrt(g*h_max)) + (h0 - h_min)*sqrt(0.5*g*(1.0/h0 + 1.0/h_min));
                slope = (F_max - F0)/(h_max - h_min);
                h0 = h0 - F0/slope;
            }
            hm = h0;
            sqrtgh2 = sqrt(g*hm);
            if (hL > hR) {
                sqrtgh1 = sqrt(g*hL);
                /* Eqn (13.55) in the FVMHP book */
                um = uL + 2.0*sqrtgh1 - 2.0*sqrtgh2;
                s1m = uL + 2.0*sqrtgh1 - 3.0*sqrtgh2;
                s2m = uL + 2.0*sqrtgh1 - sqrtgh2;

                rare1 = true;
                rare2 = false;
            } else {
                sqrtgh1 = sqrt(g*hR);
                um = uR - 2.0*sqrtgh1 + 2.0*sqrtgh2;
                s1m = uR - 2.0*sqrtgh1 + sqrtgh2;
                s2m = uR - 2.0*sqrtgh1 + 3.0*sqrtgh2;
                rare1 = false;
                rare2 = true;
            }
        }
   }

}