/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FLOOD_SPEED_USER_H
#define FLOOD_SPEED_USER_H
// #include <fc2d_cudaclaw.h>
#include <fc2d_geoclaw.h>
#include <fc2d_cudaclaw_options.h>
#include <cudaclaw_user_fort.h>

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include "riemann_source/variables.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int cuda;
    int example;
    double gravity;
    double dry_tolerance;
    double earth_radius;
    double coordinate_system;
    int mcapa;
    int is_registered;
} user_options_t;


// --- will call the riemann solvers here ----
void flood_speed_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);
void flood_speed_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2);
void flood_speed_assign_speeds(cudaclaw_cuda_speeds_t *speeds);
void flood_speed_assign_b4step2(cudaclaw_cuda_b4step2_t *b4step2);

#define GET_B4STEP2_PARAMETERS FCLAW_F77_FUNC(get_b4step2_parameters,GET_B4STEP2_PARAMETERS)
void GET_B4STEP2_PARAMETERS(const int* num_dtopo, const int* aux_finalized, double t0dtopo[], 
                            double tfdtopo[], const double* dt_max_dtopo, const double* NEEDS_TO_BE_DEFINED, 
                            const bool* variable_friction, const int* friction_index, const double* xupper, 
                            const double* yupper, const double* xlower, const double* ylower,
                            const int* test_topograpghy, const int* mtopofiles, double topowork[], double xlowtopo[], double ylowtopo[], double xhitopo[], double yhitopo[], double dxtopo[], double dytopo[], int mxtopo[], int mytopo[], int mtopoorder[], int i0topo[], int mtopo[], const int* mtoposize);

void setprob_cuda();

// --------------------------------------------
void flood_speed_link_solvers(fclaw2d_global_t *glob);
user_options_t* flood_speed_options_register (fclaw_app_t * app,
                                          const char *configfile);
void flood_speed_options_store (fclaw2d_global_t* glob, user_options_t* user);
user_options_t* flood_speed_get_options(fclaw2d_global_t* glob);

// ------------------------------------------ non-cuda functions ----------------------------
#if 0
#define FLOOD_SPEED_SETPROB FCLAW_F77_FUNC(flood_speed_setprob, FLOOD_SPEED_SETPROB)
void FLOOD_SPEED_SETPROB(const double* grav, const double* drytol,
                         const double* earth_radius, const int* coord_system,
                         const int* mcapa);
#endif
// #define FLOOD_SPEED_QINIT  FCLAW_F77_FUNC(flood_speed_qinit, FLOOD_SPEED_QINIT)

// void FLOOD_SPEED_QINIT(const int* meqn, const int* mbc,
//                     const int* mx, const int* my,
//                     const double* xlower, const double* ylower,
//                     const double* dx, const double* dy,
//                     double q[], const int* maux, double aux[]);

// //  BC (Fortran to c)
#define FLOOD_SPEED_BC2   FCLAW_F77_FUNC(flood_speed_bc2, FLOOD_SPEED_BC2)

void FLOOD_SPEED_BC2(const int* meqn, const int* mbc,
                    const int* mx, const int* my,
                    const double* xlower, const double* ylower,
                    const double* dx, const double* dy,
                    const double q[], const int* maux,
                    const double aux[], const double* t,
                    const double* dt, const int mthbc[]);


/* Mappings */
fclaw2d_map_context_t* fclaw2d_map_new_nomap();

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
