/*
Copyright (c) 2018-2022 Carsten Burstedde, Donna Calhoun, Melody Shih, Scott Aiton,
Xinsheng Qin.

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

/* GeoClaw functions */
#include "geoclaw_solver/fc2d_geoclaw.h"
#include "geoclaw_solver/fc2d_geoclaw_options.h"
#include "geoclaw_solver/fc2d_geoclaw_fort.h"
#include "geoclaw_solver/fc2d_geoclaw_output_ascii.h"
#include <fclaw_gauges.h>
#include "geoclaw_solver/fc2d_geoclaw_gauges_default.h"
#include <fclaw2d_convenience.h>  /* Needed to get search function for gauges */

/* Some mapping functions */
#include <fclaw2d_map_brick.h>
#include <fclaw2d_map.h>
#include <fclaw2d_map_query.h>

/* Needed for debugging */
#include "types.h"

#include "fc2d_cudaclaw.h"
#include "fc2d_cudaclaw_fort.h"
#include "fc2d_cudaclaw_options.h"

#include <stdlib.h>  /* For size_t */

#include <fclaw_pointer_map.h>

#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_update_single_step.h>  
#include <fclaw2d_diagnostics.h>
#include <fclaw2d_defs.h>
#include "fclaw2d_options.h"

#include <fclaw2d_patch.h>
#include <fclaw2d_clawpatch.hpp>
#include <fclaw2d_clawpatch.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_diagnostics.h>
#include <fclaw2d_clawpatch_output_ascii.h> 
#include <fclaw2d_clawpatch_output_vtk.h>
#include <fclaw2d_clawpatch_fort.h>

#include "fc2d_cudaclaw_cuda.h"  
#include "cuda_source/cudaclaw_allocate.h"   /* Needed for def. of cudaclaw_fluxes_t */

#include "fc2d_cuda_profiler.h"

// struct region_type region_type_for_debug;

/* ----------------------------- static function defs ------------------------------- */
static 
void cudaclaw_setaux(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx);

/* --------------------------- Creating/deleting patches ---------------------------- */
static
void cudaclaw_patch_setup(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    cudaclaw_setaux(glob,this_patch,this_block_idx,this_patch_idx);
}


/* --------------------- Clawpack solver functions (required) ------------------------- */

static
void cudaclaw_setprob(fclaw2d_global_t *glob)
{
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);
    if (cudaclaw_vt->fort_setprob != NULL)
    {
        cudaclaw_vt->fort_setprob();
    }
}

static
void cudaclaw_qinit(fclaw2d_global_t *glob,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx)
{
    PROFILE_CUDA_GROUP("cudaclaw_qinit",1);
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);

    FCLAW_ASSERT(cudaclaw_vt->fort_qinit != NULL); /* Must be initialized */
    int mx,my,mbc,meqn,maux,maxmx,maxmy;
    double dx,dy,xlower,ylower;
    double *q, *aux;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    // maxmx = mx;
    // maxmy = my;

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    CUDACLAW_SET_BLOCK(&this_block_idx);
    // cudaclaw_vt->fort_qinit(&maxmx,&maxmy,&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
    //                       &maux,aux);
    cudaclaw_vt->fort_qinit(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                       &maux,aux);
    CUDACLAW_UNSET_BLOCK();
}


static
void cudaclaw_bc2(fclaw2d_global_t *glob,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    double t,
                    double dt,
                    int intersects_phys_bdry[],
                    int time_interp)
{
    PROFILE_CUDA_GROUP("cudaclaw_bc2",6);
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);

    fc2d_cudaclaw_options_t *clawpack_options = fc2d_cudaclaw_get_options(glob);

    FCLAW_ASSERT(cudaclaw_vt->fort_bc2 != NULL);

    int mx,my,mbc,meqn,maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    // maxmx = mx;
    // maxmy = my;

    int *block_mthbc = clawpack_options->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[4];
    for(int i = 0; i < 4; i++)
    {
        if (intersects_phys_bdry[i])
        {
            mthbc[i] = block_mthbc[i];
        }
        else
        {
            mthbc[i] = -1;
        }
    }

    /*
      We may be imposing boundary conditions on time-interpolated data;
      and is being done just to help with fine grid interpolation.
      In this case, this boundary condition won't be used to update
      anything
    */
    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&q,&meqn);

    CUDACLAW_SET_BLOCK(&this_block_idx);
    cudaclaw_vt->fort_bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                        &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    CUDACLAW_UNSET_BLOCK();
}

/* This can be used as a value for patch_vt->patch_setup */
static
void cudaclaw_setaux(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx)
{
    PROFILE_CUDA_GROUP("cudaclaw_setaux",2);
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);

    FCLAW_ASSERT(cudaclaw_vt->fort_setaux != NULL);

    int mx,my,mbc,maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    // maxmx = mx;
    // maxmy = my;

    /* If this is a ghost patch, we only set aux values in ghost cells */
    int is_ghost = fclaw2d_patch_is_ghost(this_patch);
    int mint = 2*mbc;
    int nghost = mbc;

    CUDACLAW_SET_BLOCK(&this_block_idx);
    cudaclaw_vt->fort_setaux(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                      &maux,aux,&is_ghost,&nghost,&mint);
    CUDACLAW_UNSET_BLOCK();
}


static
void cudaclaw_b4step2(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t, double dt)

{
    PROFILE_CUDA_GROUP("cudaclaw_b4step2",1);
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);



    if (cudaclaw_vt->b4step2 != NULL)
    {
        int mx,my,mbc,meqn, maux,maxmx,maxmy;
        double xlower,ylower,dx,dy;
        double *aux,*q;
        fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
        fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

        // maxmx = mx;
        // maxmy = my;

        CUDACLAW_SET_BLOCK(&this_block_idx);
        cudaclaw_vt->fort_b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                                &dx,&dy,&t,&dt,&maux,aux);
        CUDACLAW_UNSET_BLOCK();
    }
}

static
void cudaclaw_src2(fclaw2d_global_t *glob,
                     fclaw2d_patch_t *this_patch,
                     int this_block_idx,
                     int this_patch_idx,
                     double t,
                     double dt)
{
    PROFILE_CUDA_GROUP("cudaclaw_src2",7);
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);

    int mx,my,mbc,meqn, maux,maxmx,maxmy;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    FCLAW_ASSERT(cudaclaw_vt->fort_src2 != NULL);

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    // maxmx = mx;
    // maxmy = my;

    CUDACLAW_SET_BLOCK(&this_block_idx);
    cudaclaw_vt->fort_src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                         &dx,&dy,q,&maux,aux,&t,&dt);
    CUDACLAW_UNSET_BLOCK();
}


static
double cudaclaw_update(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         void* user)
{
    PROFILE_CUDA_GROUP("cudaclaw_update",3);
    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = fc2d_cudaclaw_vt(glob);
    const fc2d_cudaclaw_options_t* cuclaw_opt;

    int iter, total, patch_buffer_len;
    size_t size, bytes;
    double maxcfl;

    /* ------------------------------- Call b4step2 ----------------------------------- */
#if 1
    if (cudaclaw_vt->b4step2 != NULL)
    {
        fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_ADVANCE_B4STEP2]);       
        cudaclaw_vt->b4step2(glob,
                           this_patch,
                           this_block_idx,
                           this_patch_idx,t,dt);
        fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_ADVANCE_B4STEP2]);       
    }
#endif

    /* -------------------------------- Main update ----------------------------------- */
    fclaw2d_timer_start_threadsafe (&glob->timers[FCLAW2D_TIMER_ADVANCE_STEP2]);  

    cuclaw_opt = fc2d_cudaclaw_get_options(glob);
    maxcfl = 0.0;


    fclaw2d_single_step_buffer_data_t *buffer_data = 
              (fclaw2d_single_step_buffer_data_t*) user;

    patch_buffer_len = cuclaw_opt->buffer_len;
    iter = buffer_data->iter;
    total = buffer_data->total_count; 
    
    /* Be sure to save current step! */
    fclaw2d_clawpatch_save_current_step(glob, this_patch);

    cudaclaw_patch_data_t* patch_data = (cudaclaw_patch_data_t*) buffer_data->user;
    maxcfl = 0;
    if (iter == 0)
    {
        /* Create array to store pointers to patch data */
        patch_data = (cudaclaw_patch_data_t*) FCLAW_ALLOC(cudaclaw_patch_data_t,1);
        size = (total < patch_buffer_len) ? total : patch_buffer_len;
        bytes = size*sizeof(cudaclaw_fluxes_t);
        
        if (cuclaw_opt->src_term > 0)
        {
            patch_data->patch_array = FCLAW_ALLOC(fclaw2d_patch_t*,size);
            patch_data->patchno_array = FCLAW_ALLOC(int,size);
            patch_data->blockno_array = FCLAW_ALLOC(int,size);
        }

        
        patch_data->flux_array = FCLAW_ALLOC(cudaclaw_fluxes_t,size); // Is it bytes or size?
        // buffer_data->user = FCLAW_ALLOC(cudaclaw_fluxes_t,bytes);
        buffer_data->user = patch_data;
    } 

    /* Buffer pointer to fluxes */
    cudaclaw_store_buffer(glob,this_patch,this_patch_idx,this_block_idx,total,iter,
                            patch_data->flux_array,
                            patch_data->patch_array,
                            patch_data->patchno_array,
                            patch_data->blockno_array);

    /* Update all patches in buffer if :
          (1) we have filled the buffer, or 
          (2) we have a partially filled buffer, but no more patches to update */

    if ((iter+1) % patch_buffer_len == 0)
    {
        /* (1) We have filled the buffer */
        maxcfl = cudaclaw_step2_batch(glob,(cudaclaw_fluxes_t*) patch_data->flux_array,
                                      patch_buffer_len,t,dt);
    }
    else if ((iter+1) == total)
    {        
        /* (2) We have a partially filled buffer, but are done with all the patches 
            that need to be updated.  */
        maxcfl = cudaclaw_step2_batch(glob,(cudaclaw_fluxes_t*) patch_data->flux_array,
                                      total%patch_buffer_len,t,dt); 
    }  

    
    /* -------------------------------- Source term ----------------------------------- */
    // Check if we have stored all the patches in the buffer
    if (((iter+1) % patch_buffer_len == 0) || ((iter+1) == total))
    {
        if (cuclaw_opt->src_term > 0)
        {   
            FCLAW_ASSERT(cudaclaw_vt->src2 != NULL);
            // iterate over patches in buffer and call src2 to update them
            for (int i = 0; i < (iter+1); i++)
            {
                cudaclaw_vt->src2(glob,
                                  patch_data->patch_array[i],
                                  patch_data->blockno_array[i],
                                  patch_data->patchno_array[i],
                                  t,dt);
            }
            FCLAW_FREE(patch_data->patch_array);
            FCLAW_FREE(patch_data->patchno_array);
            FCLAW_FREE(patch_data->blockno_array);
        }
    }   

    if (iter == total-1)
    {
        // FCLAW_FREE(patch_data->patch_array);
        FCLAW_FREE(patch_data->flux_array);   
        FCLAW_FREE(buffer_data->user);                                   
    }

    fclaw2d_timer_stop_threadsafe (&glob->timers[FCLAW2D_TIMER_ADVANCE_STEP2]);    

    return maxcfl;
}

/* ---------------------------------- Output functions -------------------------------- */

static
void cudaclaw_output(fclaw2d_global_t *glob, int iframe)
{
    const fc2d_cudaclaw_options_t* clawpack_options;
    clawpack_options = fc2d_cudaclaw_get_options(glob);

    if (clawpack_options->ascii_out != 0)
    {
        // fclaw2d_clawpatch_output_ascii(glob,iframe);
        fc2d_geoclaw_output_ascii(glob,iframe); 
    }

    if (clawpack_options->vtk_out != 0)
    {
        fclaw2d_clawpatch_output_vtk(glob,iframe);
    }

}

/* Regridding functions */
static
int cudaclaw_patch_tag4refinement(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *this_patch,
                                    int blockno,
                                    int patchno,
                                    int initflag)
{
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;
    double *q, *aux;
    int tag_patch;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    
    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
 
    int level = this_patch->level;
    double t = glob->curr_time;

    /* First check to see if we are forced to refine based on regions 
       If patch intersects a region (including time interval), this routine
       returns :  

          -- level >= maximum level allowed by any region 
             this patch intersects with. (tag_patch = 0)

          -- level < minimum level required by any region
             this patch intersects with. (tag_patch = 1)

        Otherwise, tag_patch = -1 and we should refine using usual criteria.
    */
    double xupper = xlower + mx*dx;
    double yupper = ylower + my*dy;
    int refine = 1; /* We are tagging for refinement */
    FC2D_GEOCLAW_TEST_REGIONS(&level,&xlower,&ylower,&xupper,&yupper,
                              &t,&refine, &tag_patch);

    if (tag_patch < 0)
    {
        /* Need maxlevel to get length speed_tolerance*/
        const fclaw_options_t * fclaw_opt = fclaw2d_get_options(glob);
        int maxlevel = fclaw_opt->maxlevel;
        FC2D_GEOCLAW_FORT_TAG4REFINEMENT(&mx,&my,&mbc,&meqn,&maux,&xlower,&ylower,
                                         &dx,&dy,&t,&blockno,q,aux,&level,&maxlevel,
                                         &initflag,&tag_patch);
    }

    return tag_patch;

}

static 
int cudaclaw_patch_tag4coarsening(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *fine_patches,
                                    int blockno,
                                    int patchno,
                                    int initflag)
{
    int mx,my,mbc,meqn,maux;
    double xlower[4],ylower[4],dx,dy;
    double *q[4], *aux[4];
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(glob,&fine_patches[igrid],&q[igrid],&meqn);
        fclaw2d_clawpatch_aux_data(glob,&fine_patches[igrid],&aux[igrid],&maux);

        fclaw2d_clawpatch_grid_data(glob,&fine_patches[igrid], &mx,&my,&mbc,
                                    &xlower[igrid],&ylower[igrid],&dx,&dy);
    }

    int level = fine_patches[0].level;
    double t = glob->curr_time;
    int tag_patch;

    /* Test parent quadrant : If any of the four sibling patches are in the 
       region, we consider that an intersection.  Assume Morton ordering
       on the sibling patches (0=ll, 1=lr, 2=ul, 3=ur) */
    double xupper = xlower[1] + mx*dx;
    double yupper = ylower[2] + my*dy;
    int refine = 0; /* We are tagging for coarsening */
    FC2D_GEOCLAW_TEST_REGIONS(&level,&xlower[0],&ylower[0],&xupper,&yupper,
                              &t,&refine, &tag_patch);

    if (tag_patch < 0)
    {
        /* Region tagging is inconclusive */
        const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
        int maxlevel = fclaw_opt->maxlevel;

         FC2D_GEOCLAW_FORT_TAG4COARSENING(&blockno,&mx,&my,&mbc,&meqn,&maux,xlower,ylower,
                                         &dx,&dy, &t,q[0],q[1],q[2],q[3],
                                         aux[0],aux[1],aux[2],aux[3],
                                         &level,&maxlevel, &initflag, &tag_patch);
    }

    return tag_patch;
}

static
void cudaclaw_interpolate2fine(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *coarse_patch,
                                fclaw2d_patch_t *fine_patches,
                                int blockno,
                                int coarse_patchno,
                                int fine0_patchno)
{
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;

    fclaw2d_clawpatch_grid_data(glob,coarse_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *qcoarse;
    fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    double *auxcoarse;
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_patch_t* fine_patch = &fine_patches[igrid];

        double *qfine;
        fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

        double *auxfine;
        fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

        FC2D_GEOCLAW_FORT_INTERPOLATE2FINE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                           &maux,auxcoarse,auxfine, &igrid);
    }
}

static
void cudaclaw_average2coarse(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *fine_patches,
                                fclaw2d_patch_t *coarse_patch,
                                int blockno,
                                int fine0_patchno,
                                int coarse_patchno)
{
    /* Only mx, my are needed here */
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,coarse_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *qcoarse;
    fclaw2d_clawpatch_soln_data(glob,coarse_patch,&qcoarse,&meqn);

    double *auxcoarse;
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (int igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_patch_t* fine_patch = &fine_patches[igrid];

        double *qfine;
        fclaw2d_clawpatch_soln_data(glob,fine_patch,&qfine,&meqn);

        double *auxfine;
        fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

        const fc2d_cudaclaw_options_t* clawopt = fc2d_cudaclaw_get_options(glob);
        int mcapa = clawopt->mcapa;
        FC2D_GEOCLAW_FORT_AVERAGE2COARSE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                         &maux,auxcoarse,auxfine,&mcapa,&igrid);
    }                           
}                              

/* ------------------------- Ghost filling - patch specific ------------------------ */

void cudaclaw_average_face(fclaw2d_global_t *glob,
                          fclaw2d_patch_t *coarse_patch,
                          fclaw2d_patch_t *fine_patch,
                          int idir,
                          int iface_coarse,
                          int p4est_refineFactor,
                          int refratio,
                          int time_interp,
                          int igrid,
                          fclaw2d_patch_transform_data_t* transform_data)
{
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;

    fclaw2d_clawpatch_grid_data(glob,coarse_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *qcoarse;
    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double* qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    /* These will be empty for non-manifords cases */
    double *auxcoarse;
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double *auxfine;
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    const fc2d_cudaclaw_options_t* clawopt = fc2d_cudaclaw_get_options(glob);
    int mcapa = clawopt->mcapa;

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int manifold = fclaw_opt->manifold;
    if (manifold != 0)
    {
        fclaw_global_essentialf("cpu_cudaflood : Manifold case is not handled explicit" \
                                  "in GeoFlood.");
        exit(0);
    }

    FC2D_GEOCLAW_FORT_AVERAGE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                   &maux,auxcoarse,auxfine,&mcapa,
                                   &idir,&iface_coarse,
                                   &igrid,&transform_data);
}

void cudaclaw_interpolate_face(fclaw2d_global_t *glob,
                              fclaw2d_patch_t *coarse_patch,
                              fclaw2d_patch_t *fine_patch,
                              int idir,
                              int iside,
                              int p4est_refineFactor,
                              int refratio,
                              int time_interp,
                              int igrid,
                              fclaw2d_patch_transform_data_t* transform_data)
{
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;

    fclaw2d_clawpatch_grid_data(glob,coarse_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *qcoarse;
    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    double* auxcoarse;
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double* auxfine;
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    FC2D_GEOCLAW_FORT_INTERPOLATE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,&maux,
                                       auxcoarse,auxfine, &idir, &iside,
                                       &igrid, &transform_data);
}                             


void cudaclaw_average_corner(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *coarse_patch,
                            fclaw2d_patch_t *fine_patch,
                            int coarse_blockno,
                            int fine_blockno,
                            int coarse_corner,
                            int time_interp,
                            fclaw2d_patch_transform_data_t* transform_data)
{
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,coarse_patch, &mx,&my,&mbc,
                                 &xlower,&ylower,&dx,&dy);

    double *qcoarse;
    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    double* auxcoarse;
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double* auxfine;
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    const fc2d_cudaclaw_options_t* clawopt = fc2d_cudaclaw_get_options(glob);
    int mcapa = clawopt->mcapa;

    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    int manifold = fclaw_opt->manifold;
    if (manifold != 0)
    {
        fclaw_global_essentialf("cpu_cudaflood : Manifold case is not handled explicit" \
                                  "in GeoFlood.");
        exit(0);
    }

    FC2D_GEOCLAW_FORT_AVERAGE_CORNER(&mx,&my,&mbc,&meqn,
                                     qcoarse,qfine,&maux,auxcoarse,auxfine,
                                     &mcapa,&coarse_corner,
                                     &transform_data);
}                           

void cudaclaw_interpolate_corner(fclaw2d_global_t* glob,
                                fclaw2d_patch_t* coarse_patch,
                                fclaw2d_patch_t* fine_patch,
                                int coarse_blockno,
                                int fine_blockno,
                                int coarse_corner,
                                int time_interp,
                                fclaw2d_patch_transform_data_t* transform_data)

{
    int mx,my,mbc,meqn,maux;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,coarse_patch, &mx,&my,&mbc,
                                 &xlower,&ylower,&dx,&dy);

    double *qcoarse;
    fclaw2d_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw2d_clawpatch_get_q(glob,fine_patch);

    double* auxcoarse;
    fclaw2d_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);

    double* auxfine;
    fclaw2d_clawpatch_aux_data(glob,fine_patch,&auxfine,&maux);

    FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER(&mx,&my,&mbc,&meqn,
                                         qcoarse,qfine,&maux,
                                         auxcoarse,auxfine,
                                         &coarse_corner,&transform_data);
}

/* --------------------------- Parallel ghost patches -------------------------------- */

void cudaclaw_remote_ghost_setup(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *patch,
                                int blockno,
                                int patchno)
{
    fclaw2d_clawpatch_options_t* clawpatch_options;
    clawpatch_options = fclaw2d_clawpatch_get_options(glob);

    if (clawpatch_options->ghost_patch_pack_aux)
    {
       cudaclaw_setaux(glob,patch,blockno,patchno);
    }
    else
    {
         /* the aux array data has been packed and transferred as MPI messages */
    }
}

static
void cudaclaw_local_ghost_pack_aux(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *patch,
                                  int mint,
                                  double *auxpack,
                                  int auxsize, int packmode,
                                  int* ierror)
 {
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    fclaw2d_clawpatch_grid_data(glob,patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *aux;
    fclaw2d_clawpatch_aux_data(glob,patch,&aux,&maux);
    FC2D_GEOCLAW_LOCAL_GHOST_PACK_AUX(&mx,&my,&mbc,&maux,
                                          &mint,aux,auxpack,&auxsize,
                                          &packmode,ierror);
 }                                 

/* ------------------------------ Misc access functions ----------------------------- */

/* Called from application routines */
void fc2d_cudaclaw_module_setup(fclaw2d_global_t *glob)
{
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    const fc2d_cudaclaw_options_t *clawopt = fc2d_cudaclaw_get_options(glob);

    FC2D_GEOCLAW_SET_MODULES(&clawopt->mwaves,
                                &clawopt->mcapa,
                                &clawpatch_opt->meqn,
                                &clawpatch_opt->maux,
                                clawopt->mthlim,
                                clawopt->method,
                                &fclaw_opt->ax,
                                &fclaw_opt->bx,
                                &fclaw_opt->ay,
                                &fclaw_opt->by);    
}

/* ---------------------------------- Virtual table  ---------------------------------- */

static
fc2d_cudaclaw_vtable_t* cudaclaw_vt_new()
{
    return (fc2d_cudaclaw_vtable_t*) FCLAW_ALLOC_ZERO (fc2d_cudaclaw_vtable_t, 1);
}

static
void cudaclaw_vt_destroy(void* vt)
{
    FCLAW_FREE (vt);
}

fc2d_cudaclaw_vtable_t* fc2d_cudaclaw_vt(fclaw2d_global_t *glob)
{
	fc2d_cudaclaw_vtable_t* cudaclaw_vt = (fc2d_cudaclaw_vtable_t*) 
	   							fclaw_pointer_map_get(glob->vtables, "fc2d_cudaclaw");
	FCLAW_ASSERT(cudaclaw_vt != NULL);
	FCLAW_ASSERT(cudaclaw_vt->is_set != 0);
	return cudaclaw_vt;
}


void fc2d_cudaclaw_solver_initialize(fclaw2d_global_t* glob)
{
    fclaw_options_t* fclaw_opt = fclaw2d_get_options(glob);
	fclaw2d_clawpatch_options_t* clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
	fc2d_cudaclaw_options_t* clawopt = fc2d_cudaclaw_get_options(glob);

    clawopt->method[6] = clawpatch_opt->maux;

    /* We have to do this so that we know how to size the ghost patches */
    if (clawpatch_opt->ghost_patch_pack_aux)
    {
        fclaw_opt ->ghost_patch_pack_extra = 1; /* Pack the bathymetry */
        fclaw_opt->ghost_patch_pack_numextrafields = clawpatch_opt->maux;
    }

    int claw_version = 5;
    fclaw2d_clawpatch_vtable_initialize(glob, claw_version);

    fclaw_gauges_vtable_t*  gauges_vt = fclaw_gauges_vt(glob);

    fclaw2d_vtable_t*                fclaw_vt = fclaw2d_vt(glob);
    fclaw2d_patch_vtable_t*          patch_vt = fclaw2d_patch_vt(glob);
    fclaw2d_clawpatch_vtable_t*      clawpatch_vt = fclaw2d_clawpatch_vt(glob);

    fc2d_cudaclaw_vtable_t*  cudaclaw_vt = cudaclaw_vt_new();

#if defined(_OPENMP)
    fclaw_global_essentialf("Current implementation does not allow OPENMP + CUDA\n");
    exit(0);
#endif       

    /* ForestClaw virtual tables */
    fclaw_vt->output_frame                   = cudaclaw_output;
    fclaw_vt->problem_setup                  = cudaclaw_setprob;

    /* Set basic patch operations */
    patch_vt->setup                          = cudaclaw_patch_setup;
    patch_vt->initialize                     = cudaclaw_qinit;
    patch_vt->physical_bc                    = cudaclaw_bc2;
    patch_vt->single_step_update             = cudaclaw_update; /* Includes b4step2 and src2*/

    /* Regridding */
    patch_vt->tag4refinement                 = cudaclaw_patch_tag4refinement;
    patch_vt->tag4coarsening                 = cudaclaw_patch_tag4coarsening;
    patch_vt->interpolate2fine               = cudaclaw_interpolate2fine;
    patch_vt->average2coarse                 = cudaclaw_average2coarse;

    /* Ghost filling */
    clawpatch_vt->fort_copy_face             = FC2D_GEOCLAW_FORT_COPY_FACE;
    clawpatch_vt->fort_copy_corner           = FC2D_GEOCLAW_FORT_COPY_CORNER;

    /* GeoClaw needs specialized averaging and interpolation routines */
    patch_vt->average_face                   = cudaclaw_average_face;
    patch_vt->interpolate_face               = cudaclaw_interpolate_face;
    patch_vt->average_corner                 = cudaclaw_average_corner;
    patch_vt->interpolate_corner             = cudaclaw_interpolate_corner;

    patch_vt->remote_ghost_setup             = cudaclaw_remote_ghost_setup;
    clawpatch_vt->fort_local_ghost_pack   = FC2D_GEOCLAW_LOCAL_GHOST_PACK;
    clawpatch_vt->local_ghost_pack_aux       = cudaclaw_local_ghost_pack_aux;

    /* Diagnostic functions partially implemented in clawpach */
    clawpatch_vt->fort_compute_error_norm    = FC2D_GEOCLAW_FORT_COMPUTE_ERROR_NORM;
    clawpatch_vt->fort_compute_patch_area    = FC2D_GEOCLAW_FORT_COMPUTE_PATCH_AREA;
    clawpatch_vt->fort_conservation_check    = FC2D_GEOCLAW_FORT_CONSERVATION_CHECK;
    clawpatch_vt->fort_timeinterp            = FC2D_GEOCLAW_FORT_TIMEINTERP;

    cudaclaw_vt->fort_setprob                     = NULL;
    cudaclaw_vt->fort_setaux                      = FC2D_GEOCLAW_SETAUX;
    cudaclaw_vt->fort_qinit                       = FC2D_GEOCLAW_QINIT;
    cudaclaw_vt->fort_bc2                         = FC2D_GEOCLAW_BC2;
    cudaclaw_vt->fort_b4step2                     = FC2D_GEOCLAW_B4STEP2;
    cudaclaw_vt->fort_src2                        = FC2D_GEOCLAW_SRC2;
    cudaclaw_vt->fort_rpn2                        = FC2D_GEOCLAW_RPN2;
    cudaclaw_vt->fort_rpt2                        = FC2D_GEOCLAW_RPT2;

    gauges_vt->set_gauge_data                = geoclaw_read_gauges_data_default;
    gauges_vt->create_gauge_files            = geoclaw_create_gauge_files_default;
    gauges_vt->normalize_coordinates         = geoclaw_gauge_normalize_coordinates;

    gauges_vt->update_gauge                  = geoclaw_gauge_update_default;
    gauges_vt->print_gauge_buffer            = geoclaw_print_gauges_default; 

/*-------*/

    // /* ForestClaw vtable items */
    // fclaw_vt->output_frame                   = cudaclaw_output;
    // fclaw_vt->problem_setup                  = cudaclaw_setprob;    

    // /* These could be over-written by user specific settings */
    // patch_vt->initialize                     = cudaclaw_qinit;
    // patch_vt->setup                          = cudaclaw_setaux; 
    // patch_vt->physical_bc                    = cudaclaw_bc2;
    // patch_vt->single_step_update             = cudaclaw_update;

    // /* Set user data */
    // patch_vt->create_user_data  = cudaclaw_allocate_fluxes;
    // patch_vt->destroy_user_data = cudaclaw_deallocate_fluxes;

    // /* Wrappers so that user can change argument list */
    // cudaclaw_vt->b4step2        = cudaclaw_b4step2;
    // cudaclaw_vt->src2           = cudaclaw_src2;

    // /* Required functions  - error if NULL */
    // cudaclaw_vt->fort_bc2       = CUDACLAW_BC2_DEFAULT;
    // cudaclaw_vt->fort_qinit     = NULL;
    // cudaclaw_vt->fort_rpn2      = NULL;
    // cudaclaw_vt->fort_rpt2      = NULL;

    // /* Optional functions - call only if non-NULL */
    // cudaclaw_vt->fort_setprob   = NULL;
    // cudaclaw_vt->fort_setaux    = NULL;
    // cudaclaw_vt->fort_b4step2   = NULL;
    // cudaclaw_vt->fort_src2      = NULL;

    cudaclaw_vt->is_set = 1;

	FCLAW_ASSERT(fclaw_pointer_map_get(glob->vtables,"fc2d_cudaclaw") == NULL);
	fclaw_pointer_map_insert(glob->vtables, "fc2d_cudaclaw", cudaclaw_vt, cudaclaw_vt_destroy);
}


/* ----------------------------- User access to solver functions --------------------------- */


/* These are here in case the user wants to call Clawpack routines directly */

/* This should only be called when a new fclaw2d_clawpatch_t is created. */
void fc2d_cudaclaw_set_capacity(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx)
{
    int mx,my,mbc,maux,mcapa;
    double dx,dy,xlower,ylower;
    double *aux, *area;
    fc2d_cudaclaw_options_t *clawopt;

    clawopt = fc2d_cudaclaw_get_options(glob);
    mcapa = clawopt->mcapa;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(glob,this_patch);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

    CUDACLAW_SET_CAPACITY(&mx,&my,&mbc,&dx,&dy,area,&mcapa,
                            &maux,aux);
}



/* -------------------------- Public interface to Clawpack wrappers --------------------*/

/* These are overkill;  it isn't obvious why the user would want these */
void fc2d_cudaclaw_setprob(fclaw2d_global_t *glob)
{
    cudaclaw_setprob(glob);
}

/* This can be set to cudaclaw_vt->src2 */
void fc2d_cudaclaw_src2(fclaw2d_global_t* glob,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    cudaclaw_src2(glob,this_patch,this_block_idx,this_patch_idx,t,dt);
}


void fc2d_cudaclaw_setaux(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    cudaclaw_setaux(glob,this_patch,this_block_idx,this_patch_idx);
}


void fc2d_cudaclaw_qinit(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    cudaclaw_qinit(glob,this_patch,this_block_idx,this_patch_idx);
}

void fc2d_cudaclaw_b4step2(fclaw2d_global_t* glob,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    cudaclaw_b4step2(glob,this_patch,this_block_idx,this_patch_idx,t,dt);
}

void fc2d_cudaclaw_bc2(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         int intersects_bc[],
                         int time_interp)
{
    cudaclaw_bc2(glob,this_patch,this_block_idx,this_patch_idx,t,dt,
                   intersects_bc,time_interp);
}













