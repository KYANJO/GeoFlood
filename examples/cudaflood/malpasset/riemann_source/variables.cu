#include <fc2d_cudaclaw.h>
#include <fc2d_cudaclaw_check.h>
#include <fc2d_cudaclaw_options.h>
#include <cudaclaw_user_fort.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_include_all.h>

__constant__ double s_grav;
__constant__ double dry_tolerance;
__constant__ double earth_radius;
__constant__ int coordinate_system;
__constant__ int mcapa;

void setprob_cuda()
{
    double grav;
    double drytol;
    double earth_rad;
    int coordinate_system_;
    int mcapa_;
    FILE *f = fopen("setprob.data","r");
    fscanf(f,"%lf",&grav);
    fscanf(f,"%lf",&dry_tolerance);
    fscanf(f,"%lf",&earth_rad);
    fscanf(f,"%d",&coordinate_system_);
    fscanf(f,"%d",&mcapa_);
    fclose(f);

    CHECK(cudaMemcpyToSymbol(s_grav, &grav, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(dry_tolerance, &drytol, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(earth_radius, &earth_rad, sizeof(double)));
    CHECK(cudaMemcpyToSymbol(coordinate_system, &coordinate_system_, sizeof(int)));
    CHECK(cudaMemcpyToSymbol(mcapa, &mcapa_, sizeof(int)));
}
