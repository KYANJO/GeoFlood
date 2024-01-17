#ifndef VARIABLES_H
#define VARIABLES_H

#include "../fc2d_cudaclaw_cuda.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

// #define PI 3.14159265358979323846
// #define deg2rad (PI / 180.0)
// #define rad2deg (180.0 / PI)

/* Define define and declare structures */
struct GeofloodVars{
    int mcapa;
    int coordinate_system;
    double gravity;
    double dry_tolerance;
    double earth_radius;
    double deg2rad;
};

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* VARIABLES_H */
