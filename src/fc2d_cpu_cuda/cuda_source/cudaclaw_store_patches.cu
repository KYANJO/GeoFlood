#include "../fc2d_cpucuda.h"
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fc2d_cpucuda_options.h>

#include "cudaclaw_allocate.h"
#include <fc2d_cuda_profiler.h>


void cudaclaw_store_buffer(fclaw2d_global_t* glob,
                           fclaw2d_patch_t *this_patch,
                           int this_patch_idx,
                           int this_block_idx,
                           int total, int iter,
                           cudaclaw_fluxes_t* flux_array)
{
    PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer",4);

    const fc2d_cpucuda_options_t *cuclaw_opt = fc2d_cpucuda_get_options(glob);
    fc2d_cpucuda_vtable_t*  cudaclaw_vt = fc2d_geoclaw_vt(glob);


    cudaclaw_fluxes_t *fluxes = (cudaclaw_fluxes_t*) 
               fclaw2d_patch_get_user_data(glob,this_patch);

    FCLAW_ASSERT(fluxes != NULL);

    int meqn, maux;
    double *qold_geoclaw, *aux_geoclaw;
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux_geoclaw,&maux);
    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold_geoclaw,&meqn);
   
    fluxes->qold = qold_geoclaw;
    fluxes->aux = aux_geoclaw;

    int n = iter % cuclaw_opt->buffer_len;
    flux_array[n] = *fluxes;
          
}
