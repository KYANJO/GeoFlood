#include "../flood_speed_user.h"
#include <math.h>
#include <fc2d_cudaclaw_check.h>
#include "variables.h"

/* Extern declarations*/
extern __constant__ GeofloodVars d_geofloodVars;
extern __constant__ TopoVars d_topoVars;
extern __constant__ FrictionVars d_frictionVars;
extern __constant__ AmrVars d_amrVars;

/* Function imports */
__device__ void check4nans(int meqn, int mbc, int mx, int my, double q[], 
                            double t, int ichecknan, int i, int j);

__device__ void setaux_cuda(int mbc, int mx, int my, double xlow, 
                            double ylow, double dx, double dy, int maux, 
                            double aux[], int is_ghost_in, int nghost, int mint,int i, int j);

__device__ double fc2d_geoclaw_get_dt_max_dtopo();
__device__ void fc2d_geoclaw_get_dtopo_interval(double tmin, double tmax);
__device__ bool fc2d_geoclaw_check_dtopotime(double t, double tau);
__device__ bool ghost_invalid(int i, int j, int mx, int my, int nghost, int mint);

__device__ void flood_b4step2_cuda(int mbc, int mx, int my, 
                                    int meqn, double q[],
                                    double xlower, double ylower, 
                                    double dx, double dy, 
                                    double time, double dt, int maux, 
                                    double aux[], int i, int j)
{
    /* Access the __constant__ variables in variables.h */
    double drytol = d_geofloodVars.dry_tolerance;
    double NEEDS_TO_BE_DEFINED = d_amrVars.NEEDS_TO_BE_DEFINED;

    /* check  for NaNs in the solution */
    check4nans(meqn, mbc, mx, my, q, time, 1, i, j);

    /* - check for h < 0 and reset to zero 
       - check for h < drytol
       - set hu = hv = 0 in all these cells
    */
    if (q[0] < drytol)
    {
        q[0] = fmax(q[0], 0.0);
        q[1] = 0.0;
        q[2] = 0.0;
    }

    /* Determine if time is in dtopo interval. if so, we need to update the aux array. */
    double tau;
    bool t_in_dtopo_interval = fc2d_geoclaw_check_dtopotime(time, tau);

    if (t_in_dtopo_interval) {
        /* topo arrays might have been updated by dtopo more recently than
            aux arrays were set unless atleast 1 step taken on all levels 
        */
        aux[0] = NEEDS_TO_BE_DEFINED;
        int is_ghost = 0; /* Won't be used, if is_ghost = 0 */
        int nghost = mbc;
        int mint = 2*mbc;

        setaux_cuda(mbc,mx,my,xlower,ylower,dx,dy,maux,aux,is_ghost,nghost,mint,i,j);

    }

}

__device__ cudaclaw_cuda_b4step2_t flood_speed_b4step2 = flood_b4step2_cuda;

void flood_speed_assign_b4step2(cudaclaw_cuda_b4step2_t *b4step2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(b4step2, flood_speed_b4step2, sizeof(cudaclaw_cuda_b4step2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (flood_speed_b4step2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
}

/* check4nans */
__device__ void check4nans(int meqn, int mbc, int mx, int my, double q[], double t, int ichecknan, int i, int j)
{

    for (int m = 0; m < meqn; m++)
    {
        if (isnan(q[m])) {
            printf("CUDA SOLUTION ERROR --- ABORTING CALCULATION");
            printf("At icheknan = %d", ichecknan);
            printf("    mx,my,t:%d %d %f\n", mx, my, t);
            printf("    m,i,j:%d %d %d\n", m, i, j);
            printf("    q[%d] = ", q[m]);
        }
    }
}

__device__ double fc2d_geoclaw_get_dt_max_dtopo()
{
    /* Access the __constant__ variables in variables.h */
    double dt_max_dtopo = d_topoVars.dt_max_dtopo;

    return dt_max_dtopo;
}

__device__ void fc2d_geoclaw_get_dtopo_interval(double tmin, double tmax)
{
    /* Access the __constant__ variables in variables.h */
    int num_dtopo = d_topoVars.num_dtopo;
    double *t0dtopo = d_topoVars.t0dtopo;
    double *tfdtopo = d_topoVars.tfdtopo;

    tmin = 1.0e99;
    tmax = 0.0;
    for (int i = 0; i < num_dtopo; i++)
    {
        if (t0dtopo[i] < tmin)
        {
            tmin = t0dtopo[i];
        }
        if (tfdtopo[i] >= tmax)
        {
            tmax = tfdtopo[i];
        }
    }
}

__device__ bool fc2d_geoclaw_check_dtopotime(double t, double tau)  {
    
    double tmin, tmax;

    fc2d_geoclaw_get_dtopo_interval(tmin, tmax);

    if (tmin < tmax) {
        // dtopo time is a finite interval : dtopo interval is [tmin, tmax]
        tau = (t - tmin) / (tmax - tmin);
        if (tau >= 0.0 and tau <= 1.0) {
            return true;
        } else {
            return false;
        }
    } else {
        // tmin == tmax : dtopo interval is [tmin, \infty]
        if (t >= tmax) {
            // t is in [tmin, \infty]
            tau = 1.0e99;
            return true;
        } else {
            // We haven't yet reached the dtopo interval
            tau = -1.0;
            return false;
        }
    }
}

__device__ void setaux_cuda(int mbc, int mx, int my, double xlow, 
                            double ylow, double dx, double dy, int maux, 
                            double aux[], int is_ghost_in, int nghost, int mint,int i, int j)
{
    /* # Set auxiliary arrays
        - aux[0] = Z(x,y) topography (negative below sea level)
        - if coordinate_system=2 then lat-lon coordinates on the sphere and
            aux[1] = area ratio (capacity function --set mcapa = 2)
            aux[2] = length ratio for edge
    */
    
    /* Access the __constant__ variables in variables.h */
    int coordinate_system = d_geofloodVars.coordinate_system;
    double earth_radius = d_geofloodVars.earth_radius;
    int mcapa = d_geofloodVars.mcapa;
    bool variable_friction = d_frictionVars.variable_friction;
    int friction_index = d_frictionVars.friction_index;
    double xlower = d_amrVars.xlower;
    double ylower = d_amrVars.ylower;
    double xupper = d_amrVars.xupper;
    double yupper = d_amrVars.yupper;


    bool is_ghost = (is_ghost_in != 0);

    /* Lat-long coordinate system in use, check input variables */
    if (coordinate_system == 2)
    {
        if (mcapa != 2 || maux < 3)
        {
            printf("ERROR in setaux:  for coordinate_system==2");
            printf("       need mcapa == 2 and maux >= 3\n");
            printf("       have mcapa = %d, maux = %d\n", mcapa, maux);
            return;
        }
    }

    if (coordinate_system == 1) {
        if (mcapa > 0) {
            printf("ERROR in setaux:  for coordinate_system==1");
            printf("       mcapa == 0 and maux == 1\n");
            printf("       have mcapa = %d, maux = %d\n", mcapa, maux);
            return;
        } else if (maux > 1) {
            /* should not need to set aux[1] in this case but for some reason   it blows e.g. in bowl-radial if maux>1 
            */
            if (is_ghost && ghost_invalid(i, j, mx, my, nghost, mint)) return;
            aux[1] = 1.0;
        }
    }

    /* If using a variable friction field initialize the coefficients to 0 */
    if (variable_friction)
    {
        if (is_ghost && ghost_invalid(i, j, mx, my, nghost, mint)) return;
        aux[friction_index] = 0.0;
    }

    /* Skipped analytical topography for now */

    /* Test: compute integer indices based off same corner of domaint to reduce round off discrepancies */
    int ilo = floor((xlow - xlower + 0.50*dx) / dx);
    int jlo = floor((ylow - ylower + 0.50*dy) / dy);

    /* Set topography */
    int skipcount = 0;
    double ym = ylower + (j+jlo - 1.0)*dy;
    double yp = ylower + (j+jlo)*dy;
    double y = 0.5*(ym + yp);

    double xm = xlower + (i+ilo - 1.0)*dx;
    double xp = xlower + (i+ilo)*dx;
    double x = 0.5*(xm + xp);

    printf("in setaux %d%d%f\n",i,j,aux[0]);

    if (is_ghost && ghost_invalid(i, j, mx, my, nghost, mint)) return;

    /* set lat-long cell info */
    if (coordinate_system == 2)
    {
       aux[1] = deg2rad * earth_radius*earth_radius * (sin(yp*deg2rad) - sin(ym*deg2rad)) / dy;
       aux[2] = ym * deg2rad;
    }

    /* skip setting aux[0] in ghost cell if outside physical domain since
     topo files may not cover ghost cell, and values should be extrapolated, which is done in next set of loops*/
     if (((y>yupper) || (y<ylower)) && ((x>xupper) || (x<xlower))) return;

    /* use input topography files if available */


}

__device__ bool ghost_invalid(int i, int j, int mx, int my, int nghost, int mint)
{
    bool inner = (i > mint && i < mx - mint + 1) && (j > mint && j < my - mint + 1);
    bool outer = (i < 1 - nghost) || (i > mx + nghost) || (j < 1 - nghost) || (j > my + nghost);

    return (inner || outer);
}

