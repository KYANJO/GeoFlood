/* 
@author: Brian Kyanjo
@date: 31 July 2023
@description: Solves normal Riemann problems for the 2D shallow water equations (swe) with 
topography:
            h_t + (hu)_x + (hv)_y = 0
            (hu)_t + (hu^2 + 1/2gh^2)_x + (huv)_y = -ghb_x
            (hv)_t + (huv)_x + (hv^2 + 1/2gh^2)_y = -ghb_y
where h is the height, u is the x velocity, v is the y velocity, g is the gravitational constant, and b is the topography.
@input: ql - conatins the state vector at the left edge of each cell
        qr - contains the state vector at the right edge of each cell
        
        This data is along a slice in the x-direction if idir = 1 or along a slice in the y-direction if idir = 2.

        idir - indicates the direction of the slice

@note: - The ith Riemann problem has left state qr(i-1,:) and right state ql(i,:).
       - This solver allows the user to easily select a Riemann solver in riemann_solvers.c,    this routine initializes all the variables for the swe, accounting for wet dry boundary, dry cells, wave speeds, etc.
       
@reference: David L. George
*/

#include <fc2d_cudaclaw.h>
#include <fc2d_cudaclaw_check.h>
#include "geoflood_riemann_utils.h"

__constant__ double s_grav;
__constant__ double drytol;
__constant__ double earth_radius;
__constant__ int coordinate_system;
__constant__ int mcapa;

void setprob_cuda()
{
    double grav;
    double dry_tolerance;
    double earth_rad;
    int coordinate_system_;
    int mcapa_;
    FILE *f = fopen("setprob.data","r");
    fscanf(f,"%lf",&grav);
    fscanf(f,"%lf",&dry_tolerance);
    fscanf(f,"%lf",&earth_rad);
    fscanf(f,"%d",&coordinate_system_);
    fscanf(f,"%d",&mcapa_);
    fclose(f);

    CHECK(cudaMemcpyToSymbol(s_grav, &grav, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(drytol, &dry_tolerance, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(earth_radius, &earth_rad, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(coordinate_system, &coordinate_system_, sizeof(int)));
    CHECK(cudaMemcpyToSymbol(mcapa, &mcapa_, sizeof(int)));
}

__device__ __constant__ double pi = 4.0*atan(1.0);
__device__ __constant__ double deg2rad = pi/180.0;


__device__ void cudaflood_rpn2(int idir, int meqn, int mwaves,
                            int maux, double ql[], double qr[],
                            double auxl[], double auxr[],
                            double fwave[], double s[], 
                            double amdq[], double apdq[])
{
    // local variables
    int m,i,mw,maxiter,mu,mv;
    double wall[2];
    double fw[2][2];
    double sw[2];

    double hR,hR,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL;
    double bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat;
    double s1m,s2m;
    double hstar,hstartest,hstarHLL,sLtest,sRtest;
    double tw,dxdc;

    bool rare1, rare2;

    // --- Intializing ------------------------------------
    if((qr[0] < 0.0) || (ql[0] < 0.0) )
    {
        printf("Negative input: hl,hr=\n %f %f\n",ql[0],qr[0]);
    }

    // Initialize Riemann problem for the grid interface
    for (mw=0; mw<=mwaves; ++mw)
    {
        sw[mw] = 0.0;
        fwave[mw] = 0.0;

    }

    // set normal direction
    mu = 1+idir;
    mv = 2-idir;

    // zero (small) negative values if they exsit
    if (qr[0] < 0.0)
    {
        qr[0] = 0.0;
        qr[1] = 0.0;
        qr[2] = 0.0;
    }

    if (ql[0] < 0.0)
    {
        ql[0] = 0.0;
        ql[1] = 0.0;
        ql[2] = 0.0;
    }

    //  skip problem if in a completely dry area
    if (qr[0] <= drytol && ql[0] <= drytol) continue;

    // Riemann problem variables
    hL = qr[0];
    hR = ql[0];
    huL = qr[mu];
    huR = ql[mu];
    bL = auxr[0];
    bR = auxl[0];

    hvL = qr[mv];
    hvR = ql[mv];

    // check for wet/dry boundary
    if (hR > drytol)
    {
        uR = huR/hR;
        vR = hvR/hR;
        phiR = 0.5*s_grav*(hR*hR) + (huR*huR)/hR
    }
    else
    {
        hR = 0.0;
        huR = 0.0;
        hvR = 0.0;
        uR = 0.0;
        vR = 0.0;
        phiR = 0.0;
    }

    if (hL > drytol)
    {
        uL = huL/hL;
        vL = hvL/hL;
        phiL = 0.5*s_grav*(hL*hL) + (huL*huL)/hL
    }
    else
    {
        hL = 0.0;
        huL = 0.0;
        hvL = 0.0;
        uL = 0.0;
        vL = 0.0;
        phiL = 0.0;
    }

    wall[0] = 1.0;
    wall[1] = 1.0;
    wall[2] = 1.0;
    if (hR <= drytol)
    {
        riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,rare1,rare2,0,drytol,s_grav);
        hstartest = fmax(hK,hstar);
        if (hstartest + bL < bR) //right state should become ghost values that mirror left for wall problem
        {
            wall[1] = 0.0;
            wall[2] = 0.0;
            hR = hL;
            huR = -huL;
            bR = bL;
            phiR = phiL;
            uR = -uL;
            vR = vL;
        }
        else if (hL+bL < bR)
        {
            bR = hL +bL;
        }
    }
    else if (hL <= drytol)
    {
        riemanntype(hR,hR,uR,-uR,hstar,s1m,s2m,rare1,rare2,0,drytol,s_grav);
        hstartest = fmax(hK,hstar);
        if (hstartest + bR < bL) //left state should become ghost values that mirror right for wall problem
        {
            wall[0] = 0.0;
            wall[2] = 0.0;
            hL = hR;
            huL = -huR;
            bL = bR;
            phiL = phiR;
            uL = -uR;
            vL = vR;
        }
        else if (hR+bR < bL)
        {
            bL = hR +bR;
        }
    }

    //  determine wave speeds
    sL = uL - sqrt(s_grav*hL); // 1 wave speed of left state
    sR = uR + sqrt(s_grav*hR); // 2 wave speed of right state

    uhat = (sqrt(s_grav*hL)*uL + sqrt(s_grav*hR)*uR)/(sqrt(s_grav*hL) + sqrt(s_grav*hR)); // Roe velocity
    chat = sqrt(0.5*s_grav*(hL+hR)); // Roe speed of sound
    sRoe1 = uhat - chat; // Roe wave speed 1 wave
    sRoe2 = uhat + chat; // Roe wave speed 2 wave

    sE1 = fmin(sL,sRoe1); // Einfeldt wave speed 1 wave
    sE2 = fmax(sR,sRoe2); // Einfeldt wave speed 2 wave

    // --- end initializing ------------------------------------

    // solve Riemann problem
    maxiter = 0;

    riemann_aug_JCP(maxiter,2,2,hL,hR,huL,huR,hvL,hvR,bL,bR,
                    uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,s_grav,sw,fw);

    // eliminate ghost fluxes for wall
    for (mw=0; mw<2; mw++)
    {
        sw[mw] = sw[mw]*wall[mw];
        fw[mw] = fw[mw]*wall[mw];
    }

    for (mw=0; mw<mwaves; mw++)
    {
        s[mw] = sw[mw];
        fwave[mw] = fw[mw];
        fwave[mu] = fw[mw];
        fwave[mv] = fw[mw];
    }

    // -- Capacity for Mapping from Latitude Longitude to physical space ----
    if (mcapa > 0)
    {
        if (idir == 0)
        {
            dxdc = earth_radius*deg2rad;
        }
        else
        {
            dxdc = earth_radius*cos(auxl[2])*deg2rad;
        }
        for (mw=0; mw<mwaves; mw++)
        {
            s[mw] = dxdc*s[mw];
            fwave[mw] = dxdc*fwave[mw];
        }
    }

    // --- compute fluctuations ------------------------------------
    amdq[0:2] = 0.0;
    apdq[0:2] = 0.0;
    for (mw=0; mw<mwaves; mw++)
    {
        if (s[mw] < 0.0)
        {
            amdq[0:2] = amdq[0:2] + fwave[mw];
        }
        else if (s[mw] > 0.0)
        {
            apdq[0:2] = apdq[0:2] + fwave[mw];
        }
        else
        {
            amdq[0:2] = amdq[0:2] + 0.5*fwave[mw];
            apdq[0:2] = apdq[0:2] + 0.5*fwave[mw];
        }
    }

     


}

__device__ cudaclaw_cuda_rpn2_t geoflood_rpn2 = cudaflood_rpn2;

void geoflood_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpn2, geoflood_rpn2, sizeof(cudaclaw_cuda_rpn2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (cudaflood_rpn2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}