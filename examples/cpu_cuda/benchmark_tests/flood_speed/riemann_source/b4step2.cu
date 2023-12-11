#include "../flood_speed.h"
#include <math.h>
#include <fc2d_cudaclaw_check.h>

/* Function imports */
__device__ void check4nans(int meqn, int mbc, int mx, int my, double q[], double t, int ichecknan, int i, int j);

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
    bool t_in_dtopo_interval = 

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

__device__ void fc2d_geoclaw_get_dtopo_interval(double tmin, double tmax)
{
    
}