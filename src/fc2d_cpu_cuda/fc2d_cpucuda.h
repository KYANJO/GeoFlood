/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Yu-Hsuan Shih, Scott Aiton
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

#ifndef FC2D_CPUCUDA_H
#define FC2D_CPUCUDA_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */
#include <fc2d_cudaclaw_cuda.h>  /* Needed for cuda_rpn2, cuda_rpt2 and other cuda functions */
#include <fc2d_cpucuda_fort.h>  /* Needed for virtual functions */
#include "fc2d_cpucuda_options.h"
#include "geoflood_options_user.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#define MAXCFL_CAP 1000.0

typedef struct fc2d_cpucuda_vtable fc2d_cpucuda_vtable_t;

/* Forward declarations */
struct fclaw2d_patch_transform_data;
struct fclaw2d_global;
struct fclaw2d_patch;
struct geoclaw_gauge;
struct cudaclaw_fluxes;


void fc2d_geoclaw_run(struct fclaw2d_global *glob);

/* ------------------------------------- Access functions ---------------------------------- */

void fc2d_geoclaw_module_setup(struct fclaw2d_global *glob);

void fc2d_geoclaw_output(struct fclaw2d_global *glob, int iframe);

/* ------------------------------------- Virtual table ----------------------------------- */

/**
 * @brief Initialize the geoclaw solver
 * 
 * fclaw2d_vtables_intialize should be called before this function.
 * 
 * fclaw_options, fclaw2d_clawpatch_options, and fc2d_cpucuda_options should be stored in glob.
 * fclaw_options, and, fc2d_cpucuda_options will be changed in this call.
 * 
 * @param global the global context
 */

/* Initialize the model based on the CPU (MPI) version of the solvers */
void fc2d_geoclaw_solver_initialize(struct fclaw2d_global *glob);

/* Initialize the model based on accelerated Riemann solvers */
// void fc2d_cudaclaw_solver_initialize(struct fclaw2d_global* glob);

#define CUDACLAW_SWAP_DATA FCLAW_F77_FUNC(cudaclaw_swap_data,CUDACLAW_SWAP_DATA)
void CUDACLAW_SWAP_DATA(const int* mx, const int *my, const int *mbc, const int* meqn, const int* maux,
                        double qold_geoclaw[], double qold_cudaclaw[], double aux_geoclaw[], double aux_cudaclaw[],
                        const int* geoclaw2cudaclaw);

/**
 * @brief Get the geoclaw vtable
 * 
 * @param global the global context
 * @return fc2d_cpucuda_vtable_t* the vtable
 */
fc2d_cpucuda_vtable_t* fc2d_geoclaw_vt(struct fclaw2d_global *glob);

struct fc2d_cpucuda_vtable
{
    fc2d_geoclaw_setprob_t  setprob;
    fc2d_geoclaw_bc2_t      bc2;
    fc2d_geoclaw_qinit_t    qinit;
    fc2d_geoclaw_setaux_t   setaux;
    fc2d_geoclaw_b4step2_t  b4step2;
    fc2d_geoclaw_src2_t     src2;
    fc2d_geoclaw_rpn2_t     rpn2;
    fc2d_geoclaw_rpt2_t     rpt2;
    fc2d_geoclaw_fluxfun_t  fluxfun;

    cudaclaw_cuda_rpn2_t   cuda_rpn2;
    cudaclaw_cuda_rpt2_t   cuda_rpt2;
    cudaclaw_cuda_src2_t   cuda_src2;
    cudaclaw_cuda_b4step2_t   cuda_b4step2;
    cudaclaw_cuda_speeds_t cuda_speeds;

    int is_set;
};


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_CLAWPACH5_H */
