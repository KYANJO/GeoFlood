#include "../fc2d_geoclaw.h"
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fc2d_geoclaw_options.h>

#include "cudaclaw_allocate.h"
#include <fc2d_cuda_profiler.h>

// #include "../fc2d_cudaclaw_fort.h"

// #include "data_swap.h"

void cudaclaw_store_buffer(fclaw2d_global_t* glob,
                           fclaw2d_patch_t *this_patch,
                           int this_patch_idx,
                           int this_block_idx,
                           int count, int iter, 
                           cudaclaw_fluxes_t* flux_array,
                           fclaw2d_patch_t** patch_array,
                           int* patchno_array,
                           int* blockno_array)
{
    PROFILE_CUDA_GROUP("fc2d_cudaclaw_store_buffer",4);

    const fc2d_geoclaw_options_t *cuclaw_opt = fc2d_geoclaw_get_options(glob);
    fc2d_geoclaw_vtable_t*  cudaclaw_vt = fc2d_geoclaw_vt(glob);


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
#if 1
    if (cuclaw_opt->src_term > 0 && cudaclaw_vt->src2 != NULL)
    {
        patch_array[n] = this_patch;
        patchno_array[n] = this_patch_idx;
        blockno_array[n] = this_block_idx;
    }
#endif    
        
}
