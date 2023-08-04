/* 
@author: Brian Kyanjo
@date: 31 July 2023
@description: Solves transvere Riemann problems for the 2D shallow water equations (swe) with 
topography
*/

#include "geoflood_riemann_utils.h"
 
// --remember to call these variables from geoclaw_module.f90 ---//
__constant__ double grav;
__constant__ double dry_tolerance;
__constant__ double earth_radius;
__constant__ double deg2rad;
__constant__ int coordinate_system;

__device__ double fmax(double, double)
__device__ double fmin(double x, double y)
__device__ double fabs(double x)
__device__ double pow(double x, double y)
__device__ double sqrt(double x)
__device__ double sin(double x)
__device__ double cos(double x)


__device__ void cudaflood_rpt2(int idir, int meqn, int mwaves, int maux,
                                double ql[], double qr[], double aux1[], 
                                double aux2[], double aux3[], int imp, 
                                double asdq[], double bmasdq[], double bpasdq[]) //<-- added imp, don't forget to update it in cudaclaw_flux2.cu
{
    int i,m,mw,mu,mv;

    double s[2];
    double r[2][2];
    double beta[2];
    double abs_tol;
    double hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br;
    double uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r;
    double delf1,delf2,delf3,dxdcd,dxdcu;
    double dxdcm,dxdcp,topo1,topo3,eta;

    abs_tol = dry_tolerance;

    mu = 1+idir;
    mv = 2-idir;

    hl = qr[0];
    hr = ql[0];
    hul = qr[mu];
    hur = ql[mu];
    hvl = qr[mv];
    hvr = ql[mv];

    //--- determine velocity from momentum ---//
    if (hl < abs_tol)
    {   
        hl = 0.0;
        ul = 0.0;
        vl = 0.0;
    }
    else
    {
        ul = hul/hl;
        vl = hvl/hl;
    }

    if (hr < abs_tol)
    {
        hr = 0.0;
        ur = 0.0;
        vr = 0.0;
    }
    else
    {
        ur = hur/hr;
        vr = hvr/hr;
    }
    
    for (mw=0; mw < mwaves; mw++)
    {
        s[mw] = 0.0;
        beta[mw] = 0.0;
        for (m=0; m<meqn; m++)
        {
            r[mw][m] = 0.0;
        }
    }
    dxdcp = 1.0;
    dxdcm = 1.0;

    if (hl <= abs_tol && hr <= abs_tol) continue;

    // check and see if cell that transverse waves are going in is high and dry
    // if (imp == 0)
    // {
        eta = qr[0] + aux2[0];
        topo1 = aux1[0];
        topo3 = aux3[0];
    // }
    // else
    // {
    //     eta = ql[0] + aux2[0];
    //     topo1 = aux1[0];
    //     topo3 = aux3[0];
    // }
    if (eta < fmax(topo1,topo3)) continue;

    if (coordinate_system == 1)
    {
        if (ixy == 1 )
        {
            dxdcp = (earth_radius*deg2rad);
            dxdcm = dxdcp;
        }
        else
        {
            // if (imp == 0)
            // {
                dxdcp = earth_radius*cos(aux3[2])*deg2rad;
                dxdcm = earth_radius*cos(aux1[2])*deg2rad;
            // }
            // else
            // {
            //     dxdcp = earth_radius*cos(aux3[2])*deg2rad;
            //     dxdcm = earth_radius*cos(aux1[2])*deg2rad;
            // }
        }
    }
    // ---- determine some speeds necessary for the Jacobian ----//
    vhat = (vr*sqrt(hr))/(sqrt(hr)+sqrt(hl)) + (vl*sqrt(hl))/(sqrt(hr)+sqrt(hl));
    uhat = (ur*sqrt(hr))/(sqrt(hr)+sqrt(hl)) + (ul*sqrt(hl))/(sqrt(hr)+sqrt(hl));
    hhat = (hr + hl)/2.0;

    roe1 = vhat - sqrt(grav*hhat);
    roe3 = vhat + sqrt(grav*hhat);

    s1l = vl - sqrt(grav*hl);
    s3r = vr + sqrt(grav*hr);

    s1 = fmin(roe1,s1l);
    s3 = fmax(roe3,s3r);

    s2 = 0.5*(s1+s3);

    s[0] = s1;
    s[1] = s2;
    s[2] = s3;

    // ---- determine asdq decomposition (beta) ----//
    delf1 = asdq[0];
    delf2 = asdq[mu];
    delf3 = asdq[mv];

    beta[0] = (s3*delf1/(s3-s1)) - (delf3/(s3-s1));
    beta[1] = -s2*delf1 + delf2;
    beta[2] = (delf3/(s3-s1)) - (s1*delf1/(s3-s1));

    // --- Set-up eigenvectors matrix (r) ---//
    r[0][0] = 1.0;
    r[1][0] = s2;
    r[2][0] = s1;

    r[0][1] = 0.0;
    r[1][1] = 1.0;
    r[2][1] = 0.0;

    r[0][2] = 1.0;
    r[1][2] = s2;
    r[2][2] = s3;

    // ---- Compute fluctuations ----//
    bmasdq[0] = 0.0;
    bpasdq[0] = 0.0;
    bmasdq[1] = 0.0;
    bpasdq[1] = 0.0;
    bmasdq[2] = 0.0;
    bpasdq[2] = 0.0;
    for (mw=0; mw < mwaves; mw++)
    {
        if (sw[mw] < 0.0)
        {
            bmasdq[0] = bmasdq[0] + dxdcm*s[mw]*beta[mw]*r[0][mw];
            bmasdq[mu] = bmasdq[mu] + dxdcm*s[mw]*beta[mw]*r[1][mw];
            bmasdq[mv] = bmasdq[mv] + dxdcm*s[mw]*beta[mw]*r[2][mw];
        }
        else if (s[mw] > 0.0)
        {
            bpasdq[0] = bpasdq[0] + dxdcp*s[mw]*beta[mw]*r[0][mw];
            bpasdq[mu] = bpasdq[mu] + dxdcp*s[mw]*beta[mw]*r[1][mw];
            bpasdq[mv] = bpasdq[mv] + dxdcp*s[mw]*beta[mw]*r[2][mw];
        }
    }
}

__device__ cudaclaw_cuda_rpt2_t geoflood_rpt2 = cudaflood_rpt2;

void geoflood_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(rpt2, geoflood_rpt2, sizeof(cudaclaw_cuda_rpt2_t));

    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (cudaflood_rpt2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}