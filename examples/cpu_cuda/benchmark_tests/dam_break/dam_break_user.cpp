/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "dam_break_user.h"

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fc2d_cpucuda.h>

static
void dam_break_problem_setup(fclaw2d_global_t* glob)
{
    user_options_t* user_opt = geoflood_get_options(glob);
    fclaw2d_domain_barrier (glob->domain);
    if (user_opt->cuda != 0)
    {
        setprob_cuda();
    }
}


void dam_break_link_solvers(fclaw2d_global_t *glob)
{   
    fclaw2d_vtable_t *vt = fclaw2d_vt(glob);
    vt->problem_setup = &dam_break_problem_setup;  /* Version-independent */
    /* These are set by GeoClaw for convenience, but the user
       can set these with customized functions, if desired. */
    fc2d_cpucuda_vtable_t* geoclaw_vt = fc2d_geoclaw_vt(glob);

    geoclaw_vt->qinit = &DAM_BREAK_QINIT; /* initial conditions */
}
