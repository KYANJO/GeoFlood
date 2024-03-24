/* 
    @author Brian Kyanjo briankyanjo@u.boisestate.edu
    @date 2024.03.23
    @brief b4step2 kernel function for GeoFlood: called before each call to step
     - use to set time-dependent aux arrays or perform other tasks

     - This particular kernel routine sets negative values of q[0] to zero,
       as wells as the corresponding values of q[1] and q[2] to zero.
    - This is done to ensure that the depth is always positive.
    - This should occur only becuase of round-off errors.

    - Routines for moving topography have not been implemented yet, but will be
      added in the future.

*/

#include "../fc2d_cudaclaw_cuda.h"
#include "variables.h"
#include <math.h>
#include <fc2d_geoclaw.h>
#include <fc2d_cudaclaw_check.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_include_all.h>

/* Extern declarations*/
// extern __constant__ GeofloodVars d_geofloodVars;

__device__ void cuda_flood_b4step2(double drytol, double q[])
{
    // check for NaNs in q
    // bool hasNaNs = isnan(q[0]) || isnan(q[1]) || isnan(q[2]);
    // if (hasNaNs)
    // {
    //     printf("NaNs detected in q\n");
    //     return; // Early exit if any NaNs found to avoid further computation.
    // }

    // check for q[0] < 0 and reset to zero
    // check for q[0] < drytol and set q[1] and q[2] to zero
    if (q[0] < drytol)
    {
        q[0] = fmax(q[0], 0.0);
        q[1] = 0.0;
        q[2] = 0.0;
    }

    // q[0] = fmax(q[0], 0.0); // Ensure q[0] is not negative, applies unconditionally

    // // Calculate condition once and reuse, avoiding branching
    // double condition = (q[0] < drytol);

    // // Set q[1] and q[2] to 0 if condition is true (q[0] < drytol), otherwise leave them unchanged
    // q[1] *= (1.0 - condition);
    // q[2] *= (1.0 - condition);
}

__device__ cudaclaw_cuda_b4step2_t cudaflood_b4step2 = cuda_flood_b4step2;

void cudaflood_assign_b4step2(cudaclaw_cuda_b4step2_t *b4step2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(b4step2, cudaflood_b4step2, sizeof(cudaclaw_cuda_b4step2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (cuda_flood_b4step2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}