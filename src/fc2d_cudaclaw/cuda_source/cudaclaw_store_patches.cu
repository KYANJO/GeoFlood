#include "../fc2d_cudaclaw.h"
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>
#include <fc2d_cudaclaw_options.h>

#include "cudaclaw_allocate.h"
#include <fc2d_cuda_profiler.h>

void cudaclaw_store_buffer(fclaw2d_global_t* glob,
                           fclaw2d_patch_t *this_patch,
                           int this_patch_idx,
                           int count, int iter, 
                           cudaclaw_fluxes_t* flux_array)
{
    PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer",4);
    double *qold, *aux;
    int meqn, maux;

    const fc2d_cudaclaw_options_t *cuclaw_opt = fc2d_cudaclaw_get_options(glob);

    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
               fclaw2d_patch_get_user_data(glob,this_patch);

    FCLAW_ASSERT(fluxes != NULL);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);
    
    // qold
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;
    int mbc = clawpatch_opt->mbc;
    
    double *qold_transpose = FCLAW_ALLOC(double,(mx+2*mbc)*(my+2*mbc)*meqn);
    // Write a fortran routine that transposes data form (m,i,j)-geoclaw  to (i,j,m)-cudaclaw
    
    // aux
    double *aux_transpose = FCLAW_ALLOC(double,(mx+2*mbc)*(my+2*mbc)*maux);
     // Write a fortran routine that transposes data form (m,i,j)-geoclaw  to (i,j,m)-cudaclaw

    fluxes->qold = qold_transpose;
    fluxes->aux = aux_transpose;

    flux_array[iter % cuclaw_opt->buffer_len] = *fluxes;
}
