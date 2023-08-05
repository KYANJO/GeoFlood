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

#ifndef MALPASSET_USER_H
#define MALPASSET_USER_H

#include <fc2d_cudaclaw.h>
#include <fc2d_cudaclaw_options.h>
#include <cudaclaw_user_fort.h>

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#ifdef __cplusplus
extern "C"
{
#endif

#if 0
/* Fix syntax highlighting */
#endif

typedef struct user_options
{
    int cuda;
    double gravity;
    double dry_tolerance;
    double earth_radius;
    int coordinate_system;
    int mcapa;
} user_options_t;


/* --------------------------------------- Cuda ----------------------------------------*/

  void malpasset_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);
  void malpasset_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2);
  void malpasset_assign_speeds(cudaclaw_cuda_speeds_t *speeds);

  void setprob_cuda();

void malpasset_link_solvers(fclaw2d_global_t *glob);

user_options_t* malpasset_options_register (fclaw_app_t * app,
                                       const char *configfile);

void malpasset_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* malpasset_get_options(fclaw2d_global_t* glob);

/* ----------------------------- Conservative update ---------------------------------- */

// fclaw2d_map_context_t* fclaw2d_map_new_nomap();

// #define RPN2_CONS_UPDATE FCLAW_F77_FUNC(rpn2_cons_update,RPN2_CONS_UPDATE)

// void RPN2_CONS_UPDATE(const int* meqn, const int* maux, const int* idir, const int* iface,
//                       double q[], double aux_center[], double aux_edge[], double flux[]);



/* --------------------------------------- non-Cuda ----------------------------------------*/

#if 0
#define MALPASSET_SETPROB FCLAW_F77_FUNC(malpasset_setprob, MALPASSET_SETPROB)
void MALPASSET_SETPROB();

#define USER5_SETAUX_MANIFOLD FCLAW_F77_FUNC(user5_setaux_manifold, \
                                             USER5_SETAUX_MANIFOLD)

void USER5_SETAUX_MANIFOLD(const int* mbc,
                           const int* mx, const int* my,
                           const double* xlower, const double* ylower,
                           const double* dx, const double* dy,
                           const int* maux, double aux[],
                           double xnormals[], double xtangents[],
                           double ynormals[], double ytangents[],
                           double surfnormals[],
                           double area[]);


void malpasset_patch_setup(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx);

fclaw2d_map_context_t* fclaw2d_map_new_nomap();

fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);

fclaw2d_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                  const double alpha);

#endif


#ifdef __cplusplus
}
#endif

#endif
