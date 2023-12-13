#ifndef VARIABLES_H
#define VARIABLES_H

#include "../flood_speed_user.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define PI 3.14159265358979323846
#define deg2rad (PI / 180.0)

/* Define define and declare structures */
struct GeofloodVars{
    int mcapa;
    int coordinate_system;
    double gravity;
    double dry_tolerance;
    double earth_radius;
};

struct TopoVars{
    int num_dtopo,test_topography;
    int aux_finalized, mtopofiles, mtoposize;
    double *t0dtopo, *tfdtopo;
    double dt_max_dtopo;
    double *topowork, *xlowtopo, *ylowtopo, *xhitopo, *yhitopo, *dxtopo, *dytopo;
    int *mxtopo, *mytopo, *mtopoorder, *i0topo, *mtopo;
};

struct FrictionVars{
    bool variable_friction;
    int friction_index;
};

struct AmrVars{
    double xupper;
    double yupper;
    double xlower;
    double ylower;
    double NEEDS_TO_BE_DEFINED;
};

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* VARIABLES_H */
