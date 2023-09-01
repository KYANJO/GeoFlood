/* Description: Solves normal Riemann problems for the 2D Shallow Water equations with topography:
    h_t + (hu)_x + (hv)_y = 0
    (hu)_t + (hu^2 + 1/2gh^2)_x + (huv)_y = -ghb_x
    (hv)_t + (huv)_x + (hv^2 + 1/2gh^2)_y = -ghb_y
    where h is the water depth, u and v are the horizontal and vertical velocities, and b is the topography.

On input, ql contains the state vector at the left edge of each cell
          qr contains the state vector at the right edge of each cell. 

This data is along a slice in the x-direction if ixy = 1 and in the y-direction if ixy = 2.

Note that the i^th Riemann problem has left state qr(i-1,:) 
    and right state ql(i,:).

From the basic clawpack routines, this routine is called with 
    ql = qr

    Desinged by David George, Vancouver WA, Feb. 2009
    Accelerated by 

*/