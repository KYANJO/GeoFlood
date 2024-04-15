#include "../flood_speed_user.h"
#include <math.h>
#include <fc2d_cudaclaw_check.h>
#include "variables.h"

/* ========== Extern declarations ======================== */
extern __constant__ GeofloodVars d_geofloodVars;
extern __constant__ TopoVars d_topoVars;
extern __constant__ FrictionVars d_frictionVars;
extern __constant__ AmrVars d_amrVars;

/* =========== Function prototypes ======================== */
__device__ double fc2d_geoclaw_get_dt_max_dtopo();
__device__ void fc2d_geoclaw_get_dtopo_interval(double tmin, double tmax);
__device__ bool fc2d_geoclaw_check_dtopotime(double t, double tau);
__device__ bool ghost_invalid(int i, int j, int mx, int my, int nghost, int mint);

__device__ void check4nans(int meqn, int mbc, int mx, int my, double q[], 
                            double t, int ichecknan, int i, int j);

__device__ void setaux_cuda(int mbc, int mx, int my, double xlow, 
                            double ylow, double dx, double dy, int maux, 
                            double aux[], int is_ghost_in, int nghost, int mint,int i, int j);

// #if 0
// __device__ void cellgridintegrate(double topoint, double xim, double xcell, double xip, double yjm, 
//                                 double ycell, double yjp, double xlowtopo[], double ylowtopo[], double xhitopo[], 
//                                 double yhitopo[], double dxtopo[], double dytopo[], int mxtopo[], int mytopo[], 
//                                 int mtopo[], int i0topo[], int mtopoorder[], int mtopofiles, int mtoposize, double topo[]);

// __device__ void intersection(int indicator, double area, double xintlo, double xintc, double xinthi,
//                                     double yintlo, double yintc, double yinthi, double x1lo, double x1hi,
//                                     double y1lo, double y1hi, double x2lo, double x2hi, double y2lo, double y2hi);

// __device__ double bilinearintegral_s(double xim, double xip, double yjm, double yjp, 
//                                     double x1, double x2, double y1, double y2, 
//                                     double dxx, double dyy, double z11, double z12, 
//                                     double z21, double z22); 

// __device__ double bilinearintegral(double xim, double xip, double yjm, double yjp, 
//                                     double x1, double x2, double y1, double y2, 
//                                     double dxx, double dyy, double z11, double z12, 
//                                     double z21, double z22);

// __device__ double topointegral(double xim, double xip, double yjm, double yjp, 
//                                     double xxlow, double yylow, double dxx, double dyy, 
//                                     int mxx, int myy, double zz[], int intmethod); 
// #endif
// __device__ void set_friction_field(int mx, int my, int num_ghost, int num_aux,
//                                     double xlower, double ylower, double dx, double dy, double aux[], bool is_ghost, int nghost, int mint);

/* ================== Function definitions ========================= */
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
        #if 0
        setaux_cuda(mbc,mx,my,xlower,ylower,dx,dy,maux,aux,is_ghost,nghost,mint,i,j);
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

#if 0
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
    int test_topography = d_topoVars.test_topography;
    int mtopofiles = d_topoVars.mtopofiles;
    double *topowork = d_topoVars.topowork;
    double *xlowtopo = d_topoVars.xlowtopo;
    double *ylowtopo = d_topoVars.ylowtopo;
    double *xhitopo = d_topoVars.xhitopo;
    double *yhitopo = d_topoVars.yhitopo;
    double *dxtopo = d_topoVars.dxtopo;
    double *dytopo = d_topoVars.dytopo;
    int *mxtopo = d_topoVars.mxtopo;
    int *mytopo = d_topoVars.mytopo;
    int *mtopoorder = d_topoVars.mtopoorder;
    int *i0topo = d_topoVars.i0topo;
    int *mtopo = d_topoVars.mtopo;
    int mtoposize = d_topoVars.mtoposize;

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
    if (mtopofiles > 0 && test_topography == 0)
    {
        double topo_integral = 0.0;
        cellgridintegrate(topo_integral,xm,x,xp,ym,y,yp,xlowtopo,ylowtopo,
                        xhitopo,yhitopo,dxtopo,dytopo,mxtopo,mytopo,mtopo,i0topo,
                        mtopoorder,mtopofiles,mtoposize,topowork);
        if (coordinate_system == 2)
        {
            aux[0] = topo_integral /(dx * dy * aux[1]);
        }
        else
        {
            aux[0] = topo_integral / (dx * dy);
        }
    }

    /* Copy topo to ghost cells if outside physical domain */
    int iint, jint;

    y = ylower + (j + jlo - 0.5)*dy;
    if ((y < ylower) || (y > yupper))
    {
        if (is_ghost && ghost_invalid(i, j, mx, my, nghost, mint)) return;
        x = xlower + (i + ilo - 0.5)*dx;
        iint = i + fmax(0,ceil((xlower - x)/dx)) - fmax(0,ceil((x - xupper)/dx));
        jint = j + fmax(0,ceil((ylower - y)/dy)) - fmax(0,ceil((y - yupper)/dy));
        aux[0] = aux[(iint -1) * my + (jint-1)];

    }

    x = xlower + (i + ilo - 0.5)*dx;
    if ((x < xlower) || (x > xupper))
    {
        if (is_ghost && ghost_invalid(i, j, mx, my, nghost, mint)) return;
        y = ylower + (j + jlo - 0.5)*dy;
        iint = i + fmax(0,ceil((xlower - x)/dx)) - fmax(0,ceil((x - xupper)/dx));
        jint = j + fmax(0,ceil((ylower - y)/dy)) - fmax(0,ceil((y - yupper)/dy));
        aux[0] = aux[(iint -1) * my + (jint-1)];
    }

    /* set friction coefficient  based on a set of depth levels */
    // if (friction_index > 0){

    // }

    


}

__device__ bool ghost_invalid(int i, int j, int mx, int my, int nghost, int mint)
{
    bool inner = (i > mint && i < mx - mint + 1) && (j > mint && j < my - mint + 1);
    bool outer = (i < 1 - nghost) || (i > mx + nghost) || (j < 1 - nghost) || (j > my + nghost);

    return (inner || outer);
}

__device__ void cellgridintegrate(double topoint, double xim, double xcell, double xip, double yjm, 
                                  double ycell, double yjp, double xlowtopo[], double ylowtopo[], double xhitopo[], 
                                  double yhitopo[], double dxtopo[], double dytopo[], int mxtopo[], int mytopo[], 
                                  int mtopo[], int i0topo[], int mtopoorder[], int mtopofiles, int mtoposize, double topo[])
{
    /*  cellgridintegrate integrates a unique surface, over a rectangular cell
     defined from data from multiple regular Cartesian grids
     (using the finest data available in any region)

     The rectangle has coords:
     xim <= x <= xip, yjm <= y <= yjp, with center (x,y) = (xcell, ycell)

     The intersection (with one particular grid has coords:
     xintlo <= x <= xinthi, yintlo <= y <= yinthi, with center (x,y) = (xintc, yintc)
     The _set_ version uses a recursive strategy using the formulas for
     intersections of sets.
    */
    
    int indicator;
    double area, xmlo, xmc, xmhi, ymlo, ymc, ymhi;
    double xmmlo, xmmc, xmmhi, ymmlo, ymmc, ymmhi;
    double xmmmlo, xmmmc, xmmmhi, ymmmlo, ymmmc, ymmmhi;
    double xm4lo, xm4c, xm4hi, ym4lo, ym4c, ym4hi;
    double xm5lo, xm5c, xm5hi, ym5lo, ym5c, ym5hi;
 
    /* Intialize the integral of the surface */
    topoint = 0.0;

    /* determine the type of integration needed */
    int im = 1;
    int mfid,i0;
    double cellarea;
    /* first see if the grid cell is entirely in a fine topofile */
    for (int m =1 ; m <= mtopofiles; m++)
    {
        /* look at topofiles, from fine to coarse */
        mfid = mtopoorder[m];
        i0 = i0topo[mfid];

        /* check for intersection of cell and this topofile */
        cellarea = (xip - xim)*(yjp - yjm);
        intersection(indicator,area,xmlo,xmc,xmhi,ymlo,ymc,ymhi,xim,xip,yjm,yjp,xlowtopo[mfid],xhitopo[mfid],ylowtopo[mfid],yhitopo[mfid]);

        if (indicator == 1) // cell overlaps grid
        {
            if (area == cellarea) // cell is entirely in grid
            {
                /* intergrate surface and getout of here */
                topoint += topointegral(xmlo,xmc,xmhi,ymlo,ymc,ymhi,xlowtopo[mfid],ylowtopo[mfid],dxtopo[mfid],dytopo[mfid],mxtopo[mfid],mytopo[mfid],topo[i0],im);
                return;
            } else {
                goto label222;
            }
        }
    }
    
    label222:

    /* Grid cell integral must come from more than one topofile */
    for (int m = mtopofiles - 1; m >= 0; --m) {
        /* look for topo files overlapping this grid cell.
        By decending mtopoorder from mtopoorder(mtopofiles) to
        mtopoorder(1) we are decending from the coarsest to the finest topofile.
         by the end of the nested loops, the integral of this topofile
         will be taken over regions not covered by finer topofiles
         nested loops only account for a maximum of 4 intersecting topo files
         if 5 topofiles intersect, then an error is sent.
         */
        mfid = mtopoorder[m];

        /* note that mfid indicates the topofile number or id for the "m'th" coarsest topofile
        within this do loop, we are always integrating mfid
        the nested do loops find intersections with other topofiles
        to determin the areas, but for integration of mfid only!
        */

        /* check for intersection of grid cell and this topofile */
        intersection(indicator,area,xmlo,xmc,xmhi,ymlo,ymc,ymhi,xim,xip,yjm,yjp,xlowtopo[mfid],xhitopo[mfid],ylowtopo[mfid],yhitopo[mfid]);

        int mfidint;
        if (indicator == 1) // cell overlaps grid
        {
            i0 = i0topo[mfid];
            //  integrate surface over intersection of grid and cell
            topoint += topointegral(xmlo,xmc,xmhi,ymlo,ymc,ymhi,xlowtopo[mfid],ylowtopo[mfid],dxtopo[mfid],dytopo[mfid],mxtopo[mfid],mytopo[mfid],topo[i0],im);

            // loop through grids finer than this one and subtract any integrals over intersections
            for (int mm = m - 1; mm >= 0; --mm) {
                mfidint = mtopoorder[mm];
                // check for 2nd intersection
                intersection(indicator,area,xmmlo,xmmc,xmmhi,ymmlo,ymmc,ymmhi,xmlo,xmhi,ymlo,ymhi,xlowtopo[mfidint],xhitopo[mfidint],ylowtopo[mfidint],yhitopo[mfidint]);

                if (indicator == 1) 
                {
                    // get rid of coarser integral
                    topoint -= topointegral(xmmlo,xmmc,xmmhi,ymmlo,ymmc,ymmhi,xlowtopo[mfid],ylowtopo[mfid],dxtopo[mfid],dytopo[mfid],mxtopo[mfid],mytopo[mfid],topo[i0],im);

                    /* loop through grids finer than this one add back any
                        integrals over intersections that will later get subtracted again*/
                    for (int mmm = mm - 1; mmm >= 0; --mmm) {
                        mfidint = mtopoorder[mmm];
                        // check for 3rd intersection
                        intersection(indicator,area,xmmmlo,xmmmc,xmmmhi,ymmmlo,ymmmc,ymmmhi,xmmlo,xmmhi,ymmlo,ymmhi,xlowtopo[mfidint],xhitopo[mfidint],ylowtopo[mfidint],yhitopo[mfidint]);

                        if (indicator == 1) {
                            // add back integral
                            topoint += topointegral(xmmmlo,xmmmc,xmmmhi,ymmmlo,ymmmc,ymmmhi,xlowtopo[mfid],ylowtopo[mfid],dxtopo[mfid],dytopo[mfid],mxtopo[mfid],mytopo[mfid],topo[i0],im);

                            for (int m4 = mmm -1; m4 >= 0; --m4){
                                mfidint = mtopoorder[m4];
                                // check for 4th intersection
                                intersection(indicator,area,xm4lo,xm4c,xm4hi,ym4lo,ym4c,ym4hi,xmmmlo,xmmmhi,ymmmlo,ymmmhi,xlowtopo[mfidint],xhitopo[mfidint],ylowtopo[mfidint],yhitopo[mfidint]);
                            

                                if (indicator == 1) {
                                    //  add back
                                    topoint -= topointegral(xm4lo,xm4c,xm4hi,ym4lo,ym4c,ym4hi,xlowtopo[mfid],ylowtopo[mfid],dxtopo[mfid],dytopo[mfid],mxtopo[mfid],mytopo[mfid],topo[i0],im);

                                    for (int m5 = m4 - 1; m5 >= 0; --m5) {
                                        mfidint = mtopoorder[m5];
                                        // check for 5th intersection
                                        intersection(indicator,area,xm5lo,xm5c,xm5hi,ym5lo,ym5c,ym5hi,xm4lo,xm4hi,ym4lo,ym4hi,xlowtopo[mfidint],xhitopo[mfidint],ylowtopo[mfidint],yhitopo[mfidint]);

                                        if (indicator == 1) {
                                            printf("CELLGRIDINTEGRATE:");
                                            printf(" ERROR: 5 NESTED TOPOGRIDS. MAXIMUM OF 4 ALLOWED\n");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

__device__ void intersection(int indicator, double area, double xintlo, double xintc, double xinthi,
                        double yintlo, double yintc, double yinthi, double x1lo, double x1hi,
                        double y1lo, double y1hi, double x2lo, double x2hi, double y2lo, double y2hi)
{
    /* find the intersection of two rectangles, return the intersection
     and it's area, and indicator =1
     if there is no intersection, indicator =0
    */

    /* Compute the intersection of the two rectangles */
    xintlo = fmax(x1lo, x2lo);
    xinthi = fmin(x1hi, x2hi);
    yintlo = fmax(y1lo, y2lo);
    yinthi = fmin(y1hi, y2hi);

    /* compute the ceneter of the intersection */
    xintc = 0.5*(xintlo + xinthi);
    yintc = 0.5*(yintlo + yinthi);

    if (xinthi > xintlo && yinthi > yintlo)
    {
        /* compute the area of the intersection */
        area = (xinthi - xintlo)*(yinthi - yintlo);
        indicator = 1;
    }
    else
    {
        /* no intersection */
        area = 0.0;
        indicator = 0;
    }
}

/* Billinaer integration */
__device__ double bilinearintegral(double xim, double xip, double yjm, double yjp, 
                                    double x1, double x2, double y1, double y2, 
                                    double dxx, double dyy, double z11, double z12, 
                                    double z21, double z22) 
{
    /* bilinearintegral integrates the bilinear with values z##
          over the rectangular region xim <= x <= xip, and
                                      yjm <= y <= yjp
    */

    double xlow, ylow, xhi, yhi, area, sumxi, sumeta, a, b, c, d;

    /* integrate the portion of the bilinear intersected with the rectangle cell analytically */
    /* find limits of integral */
    xlow = fmax(xim, x1);
    ylow = fmax(yjm, y1);
    xhi = fmin(xip, x2);
    yhi = fmin(yjp, y2);

    // Find the area of integration
    area = (yhi - ylow) * (xhi - xlow);
    sumxi = (xhi + xlow - 2.0 * x1) / dxx;
    sumeta = (yhi + ylow - 2.0 * y1) / dyy;

    // Find coefficients of bilinear a*xi + b*eta + c*xi*eta + d
    a = z21 - z11;
    b = z12 - z11;
    c = z22 - z21 - z12 + z11;
    d = z11;

    // Compute the integral
    return (0.5 * (a * sumxi + b * sumeta) + 0.25 * c * sumxi * sumeta + d) * area;

}

__device__ double bilinearintegral_s(double xim, double xip, double yjm, double yjp, 
                                    double x1, double x2, double y1, double y2, 
                                    double dxx, double dyy, double z11, double z12, 
                                    double z21, double z22)
{
    /*
    bilinearintegral integrates the bilinear with values z##
          over the rectangular region xim <= x <= xip, and
                                      yjm <= y <= yjp
          integration is actually done on the surface of a sphere
    */

    /* Access the __constant__ variables in variables.h */
    double Rearth = d_geofloodVars.earth_radius;
    double r2d = rad2deg; // defined in variables.h
    double d2r = deg2rad; // defined in variables.h

    double xlo, ylo, xhi, yhi, delx, dely;
    double xdiffhi, xdifflo, ydiffhi, ydifflo, xdiff2;
    double adsinint, cbsinint;
    double a, b, c, d;

    /* Integrate the portion of the bilinear intersected with the rectangle cell analytically */
    /* find limits of integral */
    xlo = fmax(xim, x1);
    xhi = fmin(xip, x1);
    ylo = fmax(yjm, y1);
    yhi = fmin(yjp, y1);
    delx = xhi - xlo;
    dely = yhi - ylo;

    /* find terms for the integration */
    xdiffhi = xhi - x1;
    xdifflo = xlo - x1;
    ydiffhi = yhi - y1;
    ydifflo = ylo - y1;
    xdiff2 = 0.5 * (xdiffhi * xdiffhi - xdifflo * xdifflo);

    cbsinint = (r2d * cos(d2r * yhi) + ydiffhi * sin(d2r * yhi)) -
               (r2d * cos(d2r * ylo) + ydifflo * sin(d2r * ylo));

    adsinint = r2d * (sin(d2r * yhi) - sin(d2r * ylo));

    /* find coefficients of bilinear a*xi + b*eta + c*xi*eta + d */
    a = (z21 - z11) / dxx;
    b = (z12 - z11) / dyy;
    c = (z22 - z21 - z12 + z11) / (dxx * dyy);
    d = z11;

    // Compute the integral
    return ((a * xdiff2 + d * delx) * adsinint + 
            r2d * (c * xdiff2 + b * delx) * cbsinint) * (Rearth * d2r) * (Rearth * d2r);
} 

/* Topo integral */
__device__ double topointegral(double xim, double xip, double yjm, double yjp, 
                                double xxlow, double yylow, double dxx, double dyy, 
                                int mxx, int myy, double zz[], int intmethod) 
{
    /*  topointegral integrates a surface over a rectangular region
     that is the intersection with a Cartesion grid (of topography data)
     the surface integrated is defined by a piecewise bilinear through the
     nodes of the Cartesian grid.
      The rectangular intersection has coords:
      xim <= x <= xip, yjm <= y <= yjp

      The Cartesian grid has coords:
      xxlow <= x <= xxhi, yylow <= y <= yyhi, with grid cell size dxx by dyy
      and mxx by myy cells.
      */

    /* Access the __constant__ variables in variables.h */
    double coordinate_system = d_geofloodVars.coordinate_system;

    /* intialize */
    double theintegral = 0.0;

    double xxhi = xxlow + (mxx - 1) * dxx;
    double yyhi = yylow + (myy - 1) * dyy;

    /* Test for small Rounding errors */
    if ((xim - xxlow) < 0.0 || (xip - xxhi) > 0.0) {
        xim = fmax(xim, xxlow);
        xip = fmin(xip, xxhi);
    }

    if ((yjm - yylow) < 0.0 || (yjp - yyhi) > 0.0) {
        yjm = fmax(yjm, yylow);
        yjp = fmin(yjp, yyhi);
    }

    double = xip - xim;
    double = yjp - yjm;

    /* Integrate piecewise bilinear over regular region */
    if (intmethod == 1){ // use bilinear method
        double djjstart = (yjm - yylow) / dyy;
        int jjstart = (int) trunc(djjstart) + 1;

        double diistart = (xim - xxlow) / dxx;
        int iistart = (int) trunc(diistart) + 1;

        double djjend = (yjp - yylow) / dyy;
        int jjend = ceil(djjend) + 1;

        double diiend = (xip - xxlow) / dxx;    
        int iiend = ceil(diiend) + 1;

        iistart = fmax(iistart, 1);
        jjstart = fmax(jjstart, 1);
        iiend = fmin(iiend, mxx);
        jjend = fmin(jjend, myy);

        for (int jj = jjstart; jj < jjend - 1; ++jj) {
            double y1 = yylow + (jj - 1) * dyy;
            double y2 = yylow + (jj) * dyy;
            int jjz1 = myy - jj + 1;
            int jjz2 = jjz1 - 1;

            for (int ii = iistart; ii < iiend - 1; ++ii) {
                double x1 = xxlow + (ii - 1) * dxx;
                double x2 = xxlow + (ii) * dxx;

                double z11 = zz[(jjz1 - 1) * mxx + (ii - 1)];
                double z12 = zz[(jjz2 - 1) * mxx + (ii - 1)];
                double z21 = zz[(jjz1 - 1) * mxx + ii];
                double z22 = zz[(jjz2 - 1) * mxx + ii];

                if (coordinate_system == 1) { // Cartesian rectangle
                    theintegral += bilinearintegral(xim, xip, yjm, yjp, x1, x2, y1, y2, dxx, dyy, z11, z12, z21, z22);
                } else if (coordinate_system == 2) { // intergrate on sphere surface
                    theintegral += bilinearintegral_s(xim, xip, yjm, yjp, x1, x2, y1, y2, dxx, dyy, z11, z12, z21, z22);
                } else {
                    printf("ERROR in topointegral:  coordinate_system must be 1 or 2\n");
                    printf("       coordinate_system = %d\n", coordinate_system);
                }
            }
        }
    } else {
        printf("ERROR in topointegral:  intmethod = 1,2 is supported \n");
        printf("       intmethod = %d\n", intmethod);
    }
    return theintegral;
}             
#endif
/* set friction field */
// __device__ void set_friction_field(int mx, int my, int num_ghost, int num_aux,
//     double xlower, double ylower, double dx, double dy,
//     double aux[], bool is_ghost, int nghost, int mint)
// {
//     /* Access the __constant__ variables in variables.h */
//     bool variable_friction = d_frictionVars.variable_friction;
//     int friction_index = d_frictionVars.friction_index;
// //     ,
// //     FrictionRegion *friction_regions, int num_friction_regions,
// //    double sea_level,
// }