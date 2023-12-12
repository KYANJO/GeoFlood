#include "../flood_speed_user.h"
#include <math.h>
#include <fc2d_cudaclaw_check.h>

/* From setprob_cuda*/
extern __constant__ double drytol; 
extern __constant__ int aux_finalized;
extern __constant__ int num_dtopo;
extern __constant__ double  dt_max_dtopo;
extern __constant__ double *t0dtopo;
extern __constant__ double *tfdtopo;
extern __constant__ double NEEDS_TO_BE_DEFINED;
extern __constant__ int coordinate_system;
extern __constant__ int mcapa;
extern __constant__ double earth_radius;
extern __constant__ double deg2rad;     
extern __constant__ bool variable_friction;
extern __constant__ int friction_index;

/* Function imports */
__device__ void check4nans(int meqn, int mbc, int mx, int my, double q[], 
                            double t, int ichecknan, int i, int j);

__device__ void setaux_cuda(int mbc, int mx, int my, double xlow, 
                            double ylow, double dx, double dy, int maux, 
                            double aux[], int is_ghost_in, int nghost, int mint);

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
    return dt_max_dtopo;
}

__device__ void fc2d_geoclaw_get_dtopo_interval(double tmin, double tmax)
{
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
                            int is_ghost, int nghost, int mint)
{
    /* # Set auxiliary arrays
        - aux[0] = Z(x,y) topography (negative below sea level)
        - if coordinate_system=2 then lat-lon coordinates on the sphere and
            aux[1] = area ratio (capacity function --set mcapa = 2)
            aux[2] = length ratio for edge
    */
    
    extern __constant__ double xupper, xlower, yupper, ylower;

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

    /* set analytical bathymetry here if requested */
    if (test_topography > 0) {
        if (is_ghost && ghost_invalid(i, j, mx, my, nghost, mint)) return;
        aux[0] = test_topo(xlow + (i - 0.5)*dx);
    }
}

__device__ bool ghost_invalid(int i, int j, int mx, int my, int nghost, int mint)
{
    bool inner = (i > mint && i < mx - mint + 1) && (j > mint && j < my - mint + 1);
    bool outer = (i < 1 - nghost) || (i > mx + nghost) || (j < 1 - nghost) || (j > my + nghost);

    return (inner || outer);
}

__device__ double test_topo(double x) {
    extern __device__ int test_topography;
    extern __device__ double topo_location, topo_left, topo_right;
    extern __device__ double topo_x0, topo_x1, topo_x2;
    extern __device__ double topo_basin_depth, topo_shelf_slope, topo_shelf_depth, topo_beach_slope;

    double topography;

    if (test_topography == 1) {
        if (x < topo_location) {
            topography = topo_left;
        } else {
            topography = topo_right;
        }
    } else if (test_topography == 2) {
        if (x < topo_x0) {
            topography = topo_basin_depth;
        } else if (x >= topo_x0 && x < topo_x1) {
            topography = topo_shelf_slope * (x - topo_x0) + topo_basin_depth;
        } else if (x >= topo_x1 && x < topo_x2) {
            topography = topo_shelf_depth;
        } else {
            topography = topo_beach_slope * (x - topo_x2) + topo_shelf_depth;
        }
    }

    return topography;
}
