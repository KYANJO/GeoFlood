/* 
    @author Brian Kyanjo briankyanjo@u.boisestate.edu
    @date 2024.03.07
    @brief src2 kernel function for GeoFlood

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
extern __constant__ GeofloodVars d_geofloodVars;


/* ================= function prototypes ================= */
__device__ double coriolis(double y);


/* ================= src2 Kernel function ================= */
// @desc: src2 kernel function for GeoFlood, contains friction and coriolis source terms
__device__ void cuda_flood_src2(int meqn, int maux, double xlower, double ylower, double dx, 
                    double dy,    double *qr, double *auxr, double t, double dt, int i, int j)
{
    double g = d_geofloodVars.gravity;
    double dry_tolerance = d_geofloodVars.dry_tolerance;
    double earth_radius = d_geofloodVars.earth_radius;
    int coordinate_system = d_geofloodVars.coordinate_system;
    double deg2rad = d_geofloodVars.deg2rad;
    bool coriolis_forcing = d_geofloodVars.coriolis_forcing;
    bool friction_forcing = d_geofloodVars.friction_forcing;
    double friction_depth = d_geofloodVars.friction_depth;
    bool variable_friction = d_geofloodVars.variable_friction;
    int num_manning = d_geofloodVars.num_manning;
    int friction_index = d_geofloodVars.friction_index;
    double manning_coefficient = d_geofloodVars.manning_coefficent;
    double manning_break = d_geofloodVars.manning_break;


    // local variables
    double depth_tolerance = 1.0e-30;
    
    // ------------- Friction source term -----------------
    if (friction_forcing)
    {
        // Extract approximate momentum
        if (qr[0] < depth_tolerance)
        {
            qr[1] = 0.0;
            qr[2] = 0.0;
        }
        else
        {
            // Apply friction source term only if in shallower water
            double coeff;
            if (qr[0] <= friction_depth)
            {
                if (!(variable_friction))
                {
                    // for (int nman = num_manning - 1; nman >= 0; nman--) 
                    {
                        // if (auxr[0] < manning_break[nman])
                        if (auxr[0] < manning_break)
                        {
                            // coeff = manning_coefficient[nman];
                            coeff = manning_coefficient;
                        }
                    }
                }
                else
                {
                    coeff = auxr[friction_index];
                }

                // Apply friction source term
                double gamma =  sqrt(pow(qr[1],2) + pow(qr[2],2))* g * pow(coeff,2) / (pow(qr[0],7/3));
                double dgamma = 1.0 + dt * gamma;
                // printf("dgamma: %f\n", dgamma);
                qr[1] = qr[1] / dgamma;
                qr[2] = qr[2] / dgamma;
            }
        }
    }
    // ------------- End of friction source term -----------------

    // ------------- Coriolis source term -----------------------
    // printf("Coriolis source term\n", coriolis_forcing);
    if (coriolis_forcing)
    {
        double y = ylower + (j - 0.5) * dy;
        double fdt = coriolis(y) * dt; // Calculate f dependent on coordinate system

        // calculate the matrix components
        double a[2][2];
        a[0][0] = 1.0 - (0.5*fdt*fdt) + pow(fdt,4)/24.0;
        a[0][1] = fdt - pow(fdt,3)/6.0;
        a[1][0] = -fdt + pow(fdt,3)/6.0;
        a[1][1] = a[0][0]; 

        // ??
        qr[1] = qr[1] * a[0][0] + qr[2] * a[0][1];
        qr[2] = qr[1] * a[1][0] + qr[2] * a[1][1];
    }
    // ------------- End of Coriolis source term -----------------
}

__device__ cudaclaw_cuda_src2_t cudaflood_src2 = cuda_flood_src2;

void cudaflood_assign_src2(cudaclaw_cuda_src2_t *src2)
{
    cudaError_t ce = cudaMemcpyFromSymbol(src2, cudaflood_src2, sizeof(cudaclaw_cuda_src2_t));
    if(ce != cudaSuccess)
    {
        fclaw_global_essentialf("ERROR (cuda_flood_src2): %s\n",cudaGetErrorString(ce));
        exit(0);
    }
}

/* Calculate the coriolis constant f 
- if coordinate system  == 1 (Cartesian) then 
    A beta-plane approximation is used and y should be in meters
- if coordinate system == 2 (spherical) then
    Grid is in lat-lon coordinates and y should be in degrees which
    is then converted to radians
*/
__device__ double coriolis(double y)
{
    int coordinate_system = d_geofloodVars.coordinate_system;
    double deg2rad = d_geofloodVars.deg2rad;
    double theta_0 = d_geofloodVars.theta_0;
    double omega = d_geofloodVars.omega;

    // Assume beta-plane approximation and y in meters
    if (coordinate_system == 1)
    {
        double theta = y / 111000.0 * deg2rad + theta_0;
        return 2.0 * omega * (sin(theta_0) + (theta - theta_0) * cos(theta_0));
    }
    else if (coordinate_system == 2)
    {
        return 2.0 * omega * sin(y*deg2rad);
    }
    else
    {
        // Unknown coordinate system, return 0.0
        return 0.0;
    }
}