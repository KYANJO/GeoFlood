#include "../flood_speed_user.h"
#include <math.h>
#include <fc2d_cudaclaw_check.h>

#include "b4step2.h"

/* From setprob_cuda*/
extern __constant__ double drytol; 
extern __constant__ int aux_finalized;
extern __constant__ int num_dtopo;
extern __constant__ double  dt_max_dtopo;
extern __constant__ double *t0dtopo;
extern __constant__ double *tfdtopo;
extern __constant__ double NEEDS_TO_BE_DEFINED;

/* Function imports */
__device__ void check4nans(int meqn, int mbc, int mx, int my, double q[], double t, int ichecknan, int i, int j);

__device__ double fc2d_geoclaw_get_dt_max_dtopo();
__device__ void fc2d_geoclaw_get_dtopo_interval(double tmin, double tmax);
__device__ bool fc2d_geoclaw_check_dtopotime(double t, double tau);

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
#if 0
        fc2d_geoclaw_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux,is_ghost,nghost,mint);
#endif
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
            printf("SOLUTION ERROR --- ABORTING CALCULATION");
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