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

#ifndef TETON_USER_H
#define TETON_USER_H
#include <fc2d_cudaclaw.h>

#include <fclaw2d_include_all.h>
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fc2d_cudaclaw.h>
#include <fc2d_cudaclaw_options.h>
#include <cudaclaw_user_fort.h>
#include <fclaw2d_clawpatch_fort.h>
#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct user_options
{
    int example;
    int cuda;
} user_options_t;


void teton_link_solvers(fclaw2d_global_t *glob);

/* --------------------------------------- Options ----------------------------------------*/
      
user_options_t* teton_options_register (fclaw_app_t * app, const char *configfile);

void teton_options_store (fclaw2d_global_t* glob, user_options_t* user);

user_options_t* teton_get_options(fclaw2d_global_t* glob);

/* --------------------------------------- Cuda ----------------------------------------*/

// void teton_assign_rpn2(cudaclaw_cuda_rpn2_t *rpn2);
// void teton_assign_rpt2(cudaclaw_cuda_rpt2_t *rpt2);

#define TETON_QINIT   FCLAW_F77_FUNC(teton_qinit,TETON_QINIT)
void TETON_QINIT(const int* meqn,const int* mbc,
                 const int* mx, const int* my,
                 const double* xlower, const double* ylower,
                 const double* dx, const double* dy,
                 double q[], const int* maux, double aux[]);



/* Mappings */
fclaw2d_map_context_t* fclaw2d_map_new_nomap();

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
