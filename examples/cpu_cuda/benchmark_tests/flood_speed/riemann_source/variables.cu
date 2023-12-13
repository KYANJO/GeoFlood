#include "../flood_speed_user.h"
#include "variables.h"
#include <math.h>
#include <fc2d_cudaclaw_check.h>

/* Declare constant memory variables */
__constant__ GeofloodVars d_geofloodVars;
__constant__ TopoVars d_topoVars;
__constant__ FrictionVars d_frictionVars;
__constant__ AmrVars d_amrVars;

void setprob_cuda(){
    int i = 0;
    char * line = NULL, *p = NULL, *eptr;
    size_t len = 0;
    ssize_t read;
    double arr[5];
    FILE *f = fopen("setprob.data","r");

    while ((read = getline(&line, &len, f)) != -1) 
    {
        p =strtok(line, " "); // get first word
        arr[i] = strtod(p,&eptr);  // convert to double
        i++; 
    }
    fclose(f);
    free(line);

    /* === declare variables === */
    /* topo variables */
    int num_dtopo_, aux_finalized_,test_topography_,mtopofiles_,mtoposize_;
    double  dt_max_dtopo_;
    double *t0dtopo_, *tfdtopo_;
    double *topowork_, *xlowtopo_, *ylowtopo_, *xhitopo_, *yhitopo_, *dxtopo_, *dytopo_;
    int *mxtopo_, *mytopo_, *mtopoorder_, *i0topo_, *mtopo_;
    
    /* friction variables */
    int friction_index_;
    bool variable_friction_;

    /* amr variables */
    double xupper_, yupper_, xlower_, ylower_;
    double NEEDS_TO_BE_DEFINED_;

    GET_B4STEP2_PARAMETERS(&num_dtopo_, &aux_finalized_, t0dtopo_, tfdtopo_,  
                            &dt_max_dtopo_, &NEEDS_TO_BE_DEFINED_,&variable_friction_,
                            &friction_index_, &xupper_, &yupper_, &xlower_, &ylower_,
                            &test_topography_, &mtopofiles_, topowork_, xlowtopo_, 
                            ylowtopo_, xhitopo_, yhitopo_, dxtopo_, dytopo_, mxtopo_,
                            mytopo_, mtopoorder_, i0topo_, mtopo_, &mtoposize_);
    

    /* === Create and populate structures on the host === */
    GeofloodVars geofloodVars;
    geofloodVars.gravity = arr[0];
    geofloodVars.dry_tolerance = arr[1];
    geofloodVars.earth_radius = arr[2];
    geofloodVars.coordinate_system = (int) arr[3];
    geofloodVars.mcapa = (int) arr[4];

    TopoVars topoVars;
    topoVars.num_dtopo = num_dtopo_;
    topoVars.aux_finalized = aux_finalized_;
    topoVars.dt_max_dtopo = dt_max_dtopo_;
    topoVars.t0dtopo = t0dtopo_;
    topoVars.tfdtopo = tfdtopo_;
    topoVars.test_topography = test_topography_;
    topoVars.mtopofiles = mtopofiles_;
    topoVars.topowork = topowork_;
    topoVars.xlowtopo = xlowtopo_;
    topoVars.ylowtopo = ylowtopo_;
    topoVars.xhitopo = xhitopo_;
    topoVars.yhitopo = yhitopo_;
    topoVars.dxtopo = dxtopo_;
    topoVars.dytopo = dytopo_;
    topoVars.mxtopo = mxtopo_;
    topoVars.mytopo = mytopo_;
    topoVars.mtopoorder = mtopoorder_;
    topoVars.i0topo = i0topo_;
    topoVars.mtopo = mtopo_;
    topoVars.mtoposize = mtoposize_;

    FrictionVars frictionVars;
    frictionVars.variable_friction = variable_friction_;
    frictionVars.friction_index = friction_index_;

    AmrVars amrVars;
    amrVars.xupper = xupper_;
    amrVars.yupper = yupper_;
    amrVars.xlower = xlower_;
    amrVars.ylower = ylower_;
    amrVars.NEEDS_TO_BE_DEFINED = NEEDS_TO_BE_DEFINED_;

    /* === Copy structures to device (constant memory) === */
    CHECK(cudaMemcpyToSymbol(d_geofloodVars, &geofloodVars, sizeof(GeofloodVars)));
    CHECK(cudaMemcpyToSymbol(d_topoVars, &topoVars, sizeof(TopoVars)));
    CHECK(cudaMemcpyToSymbol(d_frictionVars, &frictionVars, sizeof(FrictionVars)));
    CHECK(cudaMemcpyToSymbol(d_amrVars, &amrVars, sizeof(AmrVars)));
}
