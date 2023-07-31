/* 
@author: Brian Kyanjo
@date: 31 July 2023
@description: Solves normal Riemann problems for the @D shallow water equations (swe) with 
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

__constant__ double s_grav;




__device__ void rpn2shallow(int idir, int meqn, int mwaves,
                            int maux, double ql[], double qr[],
                            double auxl[], double auxr[],
                            double wave[], double s[], 
                            double amdq[], double apdq[])
{
    // set normal direction
    int mu = 1+idir;
    int mv = 2-idir;


}