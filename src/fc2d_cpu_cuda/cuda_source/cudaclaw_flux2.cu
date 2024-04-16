/*
  Copyright (c) 2018 Carsten Burstedde, Donna Calhoun, Melody Shih, Scott Aiton, 
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


#include "../fc2d_cpucuda.h"
#include "../fc2d_cudaclaw_cuda.h"

#include "../fc2d_cudaclaw_check.h"

#include <cub/block/block_reduce.cuh>  

#include "cudaclaw_allocate.h"  /* Needed to for definition of 'fluxes' */
#include "cudaclaw_flux2.h"

#include <thrust/device_vector.h>

__constant__ int order[2];
__constant__ int mthlim[FC2D_CUDACLAW_MWAVES];
__constant__ int use_fwaves;

extern __constant__ GeofloodVars d_geofloodVars;

extern "C"
{
int cudaclaw_check_parameters(int mwaves)
{
    return mwaves <= FC2D_CUDACLAW_MWAVES;
}

void cudaclaw_set_method_parameters(int *order_in, int *mthlim_in, int mwaves, 
                                    int use_fwaves_in)
{
    CHECK(cudaMemcpyToSymbol(order,order_in,2*sizeof(int)));
    CHECK(cudaMemcpyToSymbol(use_fwaves,&use_fwaves_in,sizeof(int)));
    CHECK(cudaMemcpyToSymbol(mthlim,mthlim_in,mwaves*sizeof(int)));
}

}

/* Include this here so we don't include device code in fc2d_cudaclaw_cuda.h */
__device__ double cudaclaw_limiter(int lim_choice, double r);

static
__device__
void cudaclaw_flux2_and_update(const int mx,   const int my, 
                               const int meqn, const int mbc,
                               const int maux, const int mwaves, 
                               const int mwork, int src_term,
                               const int mcapa, const double dry_tol,
                               const double xlower, const double ylower, 
                               const double dx,     const double dy,
                               double *const qold,       double *const aux, 
                               double *const fm,         double *const fp, 
                               double *const gm,         double *const gp,
                               double *const amdq_trans, double *const apdq_trans, 
                               double *const bmdq_trans, double *const bpdq_trans,
                               double *const waves,      double *const speeds,
                               double *const dtdx1d,     double *const dtdy1d,
                               double *const maxcflblocks,
                               cudaclaw_cuda_rpn2_t rpn2,
                               cudaclaw_cuda_rpt2_t rpt2,
                               cudaclaw_cuda_b4step2_t b4step2,
                               cudaclaw_cuda_src2_t src2,
                                double t,double dt)
{
    typedef cub::BlockReduce<double,FC2D_CUDACLAW_BLOCK_SIZE> BlockReduce;
    
    __shared__ typename BlockReduce::TempStorage temp_storage;

    extern __shared__ double shared_mem[];

    double* start  = shared_mem + mwork*threadIdx.x;

    /* --------------------------------- Start code ----------------------------------- */

    __shared__ double dtdx, dtdy;
    __shared__ int xs, ys, zs;
    __shared__ int ifaces_x, ifaces_y, num_ifaces;

    double mcapa_flag = (mcapa > 0);

    if (threadIdx.x == 0)
    {
        // Initialize variables related to time and grid spacing
        dtdx = dt/dx;
        dtdy = dt/dy;

        // Compute strides
        xs = 1;
        ys = (2*mbc + mx)*xs;
        zs = (2*mbc + my)*xs*ys;

        ifaces_x = mx + 2*mbc;
        ifaces_y = my + 2*mbc;
        num_ifaces = ifaces_x*ifaces_y;
    }

    // Synchronize to ensure all threads see the initialized values
    __syncthreads();

    /* ---------------------------- b4step2 -------------------------------- */
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I = iy*ys + ix;  /* Start at lower left */
        double *const qr    = start;          /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qr[mq] = qold[I_q];  
        }         

        // if (qr[0] < dry_tol)
        // {
        //     qr[0] = fmax(qr[0], 0.0);
        //     qr[1] = 0.0;
        //     qr[2] = 0.0;
        // }

        // b4step2(dry_tol, qr);
        
        qr[0] = fmax(qr[0], 0.0); // Ensure q[0] is not negative, applies unconditionally

        // Calculate condition once and reuse, avoiding branching
        double condition = (qr[0] < dry_tol);

        qr[0] = fmax(qr[0], 0.0);

        // Set q[1] and q[2] to 0 if condition is true (q[0] < drytol), otherwise leave them unchanged
        qr[1] *= (1.0 - condition);
        qr[2] *= (1.0 - condition);


        for(int mq = 0; mq < meqn; mq++)
        {
            /* copy back qr to global memory */
            int I_q = I + mq*zs;
            qold[I_q] = qr[mq];
        }
    }      
    __syncthreads(); /* Needed to be sure all qold variables are available below */ 

    if (mcapa > 0)
    {
        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        {
        
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            int I = iy*ys + ix;  /* Start at lower left */
            int I_capa = I + (mcapa-1)*zs; 

            // dtdx1d[thread_index] = dtdx/aux[I_capa];
            dtdx1d[I] = dtdx/aux[I_capa];
            dtdy1d[I] = dtdy/aux[I_capa];

        }
        __syncthreads();
    }



    /* ---------------------------- X-sweeps -------------------------------- */
    /* x-sweep : Mimic limits set by step2 and the CFL calculation in flux2 */
    if (threadIdx.x == 0) 
    {
        ifaces_x = mx + 2*mbc - 1; 
        ifaces_y = mx + 2*mbc - 2; 
        num_ifaces = ifaces_x*ifaces_y;
    }
    __syncthreads();

    // extern __shared__ double dtdx1d[];
    // if (mcapa > 0)
    // {
    //     for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    //     {
    //         /* 
    //         0 <= ix < ifaces_x
    //         0 <= iy < ifaces_y

    //         ---> 0 <= ix < mx + 2*mbc - 1
    //         ---> 0 <= iy < mx + 2*mbc - 2
    //         ---> ys+1 <= I < (mx+2*mbc-1)*ys + (mx+2*mbc-2)
    //         ---> ys+1 <= I <= (mx+2*mbc-2)*ys + (mx+2*mbc-3)

    //         Example : mbc=2; mx=8; (ix=0,iy=0) --> I = 13  (i=0,j=0)
    //                             (ix=mx+2,iy=mx+1) --> I=10*12 + 9 = 129
    //                             (i=mx+2,j=my+1)
    //         */

    //         int ix = thread_index % ifaces_x;
    //         int iy = thread_index/ifaces_x;

    //         int iadd = mbc-1;  // Shift from corner by 1 in each direction
    //         int I = (iy + iadd)*ys + (ix + iadd); 
    //         int I_capa = I + (mcapa-1)*zs; 

    //         // dtdx1d[thread_index] = dtdx/aux[I_capa];
    //         dtdx1d[I-1] = dtdx/aux[I_capa];

    //     }
    //     __syncthreads();
    // }

    double maxcfl = 0;
    // double dtdx0, dtdy0;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {

        /* 
        0 <= ix < ifaces_x
        0 <= iy < ifaces_y

        ---> 0 <= ix < mx + 2*mbc - 1
        ---> 0 <= iy < mx + 2*mbc - 2
        ---> ys+1 <= I < (mx+2*mbc-1)*ys + (mx+2*mbc-2)
        ---> ys+1 <= I <= (mx+2*mbc-2)*ys + (mx+2*mbc-3)

        Example : mbc=2; mx=8; (ix=0,iy=0) --> I = 13  (i=0,j=0)
                               (ix=mx+2,iy=mx+1) --> I=10*12 + 9 = 129
                               (i=mx+2,j=my+1)
        */

        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;  // Shift from corner by 1 in each direction
        int I = (iy + iadd)*ys + (ix + iadd); 

        double *const qr     = start;                 /* meqn        */
        double *const auxr   = qr      + meqn;         /* maux        */
        double *const ql     = auxr   + maux;         /* meqn        */
        double *const auxl   = ql     + meqn;         /* maux        */
        double *const s      = auxl   + maux;         /* mwaves      */
        double *const wave   = s      + mwaves;       /* meqn*mwaves */
        double *const amdq   = wave   + meqn*mwaves;  /* meqn        */
        double *const apdq   = amdq   + meqn;         /* meqn        */

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qr[mq] = qold[I_q];                       /* Right */
            ql[mq] = qold[I_q - 1];    /* Left  */
        }

        for(int m = 0; m < maux; m++)
        {
            int I_aux = I + m*zs;
            auxr[m] = aux[I_aux];
            auxl[m] = aux[I_aux - 1];
        }               

        
        rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq, dry_tol, mcapa);

        // double maxfl_update_x = (ix > 0 && ix <= mx + 1) ? 1.0 : 0.0;

        // if ((ix > 0 && ix <= mx + 1 ) && (iy >= 0 && iy <= mx + 1))
        {
            for (int mq = 0; mq < meqn; mq++) 
            {
                int I_q = I + mq*zs;
                fm[I_q] = amdq[mq];
                fp[I_q] = -apdq[mq]; 
                if (order[1] > 0)
                {
                    amdq_trans[I_q] = amdq[mq];                                        
                    apdq_trans[I_q] = apdq[mq];  
                }
            }
            
            int I_capa = I + (mcapa-1)*zs; // mcapa is set to 2 for latlon cordinates (-1 due to the switch between fortran and C)
            // double dtdx_ = mcapa_flag * dtdx/aux[I_capa] + (1.0 - mcapa_flag) * dtdx;
            double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
            // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;

            for(int mw = 0; mw < mwaves; mw++)
            {
                // if ((fabs(s[mw]*dtdx_) > 5.6 )&& (fabs(s[mw]*dtdx_) < 5.63))
                // {
                //     printf("ix = %d, iy = %d, maxcfl_0 = %f, s = %f\n",ix,iy,fabs(s[mw]*dtdx_),s[mw]);
                // }

                maxcfl = fmax(maxcfl,fabs(s[mw])*dtdx_);
                // maxcfl = mcapa_flag*fmax(maxcfl,-s[mw]*dtdx0);
                // maxcfl = max(maxcfl, maxfl_update_x * fabs(s[mw] * dtdx_));

                if (order[0] == 2)
                {                    
                    int I_speeds = I + mw*zs;
                    speeds[I_speeds] = s[mw];
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        int k = mw*meqn + mq;
                        int I_waves = I + k*zs;
                        waves[I_waves] = wave[k];
                    }
                }
            }
            // dtdx0 = dtdx_; /* update dtdx for next step */
        }
        // aggregate = (iy <= my + 1) ? fmax(aggregate, maxcfl) : aggregate;

    }
    __syncthreads();
    
    if (threadIdx.x == 0)
    {
        ifaces_x = mx + 2*mbc - 2; 
        ifaces_y = mx + 2*mbc - 1; 
        num_ifaces = ifaces_x*ifaces_y;
    }

    // Need to make sure thread 0 has set the values above before continuing
    __syncthreads();

    // extern __shared__ double dtdy1d[];
    // if (mcapa > 0)
    // {
    //     for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    //     {
    //         /* 
    //         0 <= ix < ifaces_x
    //         0 <= iy < ifaces_y

    //         ---> 0 <= ix < mx + 2*mbc - 2
    //         ---> 0 <= iy < mx + 2*mbc - 1
    //         ---> ys+1 <= I <  (mx+2*mbc-1)*ys + (mx+2*mbc-2)
    //         ---> ys+1 <= I <= (mx+2*mbc-2)*ys + (mx+2*mbc-3)

    //         Example : mbc=2; mx=8; (ix=0,iy=0) --> I = 13  (i=0,j=0)
    //                             (ix=mx+1,iy=mx+2) --> I = 142
    //                             (i=mx+1,j=my+2)
    //         */

    //         int ix = thread_index % ifaces_x;
    //         int iy = thread_index/ifaces_x;

    //         int iadd = mbc-1;  // Shift from corner by 1 in each direction
    //         int I = (iy + iadd)*ys + (ix + iadd);  /* Start one cell from left/bottom edge */
    //         int I_capa = I + (mcapa-1)*zs; 

    //         // dtdy1d[thread_index] = dtdy/aux[I_capa];
    //         dtdy1d[I-ys] = dtdy/aux[I_capa];
    //     }
    //     __syncthreads();
    // }

    /* ---------------------------- Y-sweeps -------------------------------- */
    // double dtdy0 = 
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {
        /* 
        0 <= ix < ifaces_x
        0 <= iy < ifaces_y

        ---> 0 <= ix < mx + 2*mbc - 2
        ---> 0 <= iy < mx + 2*mbc - 1
        ---> ys+1 <= I <  (mx+2*mbc-1)*ys + (mx+2*mbc-2)
        ---> ys+1 <= I <= (mx+2*mbc-2)*ys + (mx+2*mbc-3)

        Example : mbc=2; mx=8; (ix=0,iy=0) --> I = 13  (i=0,j=0)
                               (ix=mx+1,iy=mx+2) --> I = 142
                               (i=mx+1,j=my+2)
        */

        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;  // Shift from corner by 1 in each direction
        int I = (iy + iadd)*ys + (ix + iadd);  /* Start one cell from left/bottom edge */

        double *const qr     = start;                 /* meqn        */
        double *const auxr   = qr     + meqn;         /* maux        */
        double *const qd     = auxr   + maux;         /* meqn        */
        double *const auxd   = qd     + meqn;         /* maux        */
        double *const s      = auxd   + maux;         /* mwaves      */
        double *const wave   = s      + mwaves;       /* meqn*mwaves */
        double *const bmdq   = wave   + meqn*mwaves;  /* meqn        */
        double *const bpdq   = bmdq   + meqn;         /* meqn        */

        /* ------------------------ Normal solve in Y direction ------------------- */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qr[mq] = qold[I_q];        /* Right */
            qd[mq] = qold[I_q - ys];   /* Down  */  
        }

        for(int m = 0; m < maux; m++)
        {
            int I_aux = I + m*zs;
            auxr[m] = aux[I_aux];
            auxd[m] = aux[I_aux - ys];
        }               


        rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, bmdq, bpdq, dry_tol, mcapa);

        // double maxfl_update_y = (iy > 0 && iy <= my+1) ? 1.0 : 0.0;

        /* Set value at bottom interface of cell I */
        // if ((iy > 0 && iy <= my+1) && (ix >= 0 && ix <= mx + 1))
        {
            for (int mq = 0; mq < meqn; mq++) 
            {
                int I_q = I + mq*zs;
                gm[I_q] = bmdq[mq];
                gp[I_q] = -bpdq[mq]; 
                if (order[1] > 0)
                {
                    bmdq_trans[I_q] = bmdq[mq];                                                   
                    bpdq_trans[I_q] = bpdq[mq];
                }
            }

            int I_capa = I + (mcapa-1)*zs; // mcapa is set to 2 for latlon cordinates (-1 due to the switch between fortran and C)
            // double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;
            double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;

            

            for(int mw = 0; mw < mwaves; mw++)
            {
                // if ((fabs(s[mw]*dtdy_) > 5.6 )&& (fabs(s[mw]*dtdy_) < 5.63))
                // {
                //     printf("ix = %d, iy = %d, maxcfl_1 = %f, s = %f\n",ix,iy,fabs(s[mw]*dtdy_),s[mw]);
                // }

                maxcfl = fmax(maxcfl,fabs(s[mw])*dtdy_);
                // maxcfl = fmax(maxcfl,(s[mw])*dtdy1d[I]);
                // maxcfl = max(maxcfl, maxfl_update_y * fabs(s[mw] * dtdy_));
                if (order[0] == 2)
                {                    
                    int I_speeds = I + (mwaves + mw)*zs;
                    speeds[I_speeds] = s[mw];
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        int I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                        waves[I_waves] = wave[mw*meqn + mq];
                    }
                }
            }
            // dtdy0 = dtdy_; /* update dtdy for next step */
        }
        // aggregate = (ix <= mx + 1) ? fmax(aggregate, maxcfl) : aggregate;
    }


    // Make sure all threads are done before computing CFL
    __syncthreads();

    double aggregate = BlockReduce(temp_storage).Reduce(maxcfl,cub::Max());
    // aggregate = BlockReduce(temp_storage).Reduce(aggregate,cub::Max());
    if (threadIdx.x == 0)
    {
        maxcflblocks[blockIdx.z] = aggregate;
    }

    __syncthreads();

/* ---------------------- Second order corrections and limiters --------------------*/  

    if (order[0] == 2)
    {
        if (threadIdx.x == 0)
        {
            ifaces_x = mx + 2*mbc - 1; 
            ifaces_y = mx + 2*mbc - 2; 
            num_ifaces = ifaces_x*ifaces_y;
        }
        __syncthreads();

        
        // dtdx0 = mcapa_flag * dtdx/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
        //         + (1.0 - mcapa_flag) * dtdx;
        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            int iadd = mbc-1;  // Shift from corner by 1 in each direction
            int I = (iy + iadd)*ys + (ix + iadd);

            // int I_capa = I + (mcapa-1)*zs; 
            // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
            // double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
            // double dtdx_ave = 0.5*(dtdx0 + dtdx_);

            // double dtdx_ave = 0.5*(dtdx1d[I-1] + dtdx1d[I]);
            double dtdx_ave = mcapa_flag * 0.5*(dtdx1d[I-1] + dtdx1d[I]) + (1.0 - mcapa_flag) * dtdx;
            
            /* ------------------------------- X-directions --------------------------- */

            double *const s = start;           /* mwaves */
            double *const wave = s + mwaves;   /* meqn*mwaves */
            for(int mw = 0; mw < mwaves; mw++)
            {                
                int I_speeds = I + mw*zs;
                s[mw] = speeds[I_speeds];

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_waves = I + (mw*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                }                        

                if (mthlim[mw] > 0)
                {
                    double wnorm2 = 0, dotl=0, dotr = 0;
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        wnorm2 += pow(wave[mq],2);

                        int I_waves = I + (mw*meqn + mq)*zs;
                        dotl += wave[mq]*waves[I_waves-1];
                        dotr += wave[mq]*waves[I_waves+1];
                    }
                    wnorm2 = (wnorm2 == 0) ? 1e-15 : wnorm2;
                                            
                    double r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;
                    double wlimitr = cudaclaw_limiter(mthlim[mw],r);  

                    for (int mq = 0; mq < meqn; mq++)
                    {
                        wave[mq] *= wlimitr;
                    }                    
                }

                for(int mq = 0; mq < meqn; mq++)
                {
                    double cqxx = (1.0 - fabs(s[mw])*dtdx_ave)*wave[mq];
                    cqxx *= (use_fwaves) ? copysign(1.,s[mw]) : fabs(s[mw]);

                    int I_q = I + mq*zs;
                    fm[I_q] += 0.5*cqxx;   
                    fp[I_q] += 0.5*cqxx;  
                    if (order[1] > 0)
                    {                         
                        amdq_trans[I_q] += cqxx;   
                        apdq_trans[I_q] -= cqxx;      /* Subtract cqxx later */                         
                    }  
                }
            }
            // dtdx0 = dtdx_; /* update dtdx for next step */
        }
    }

    __syncthreads();

    if (order[0] == 2)
    {
        if (threadIdx.x == 0){
            ifaces_x = mx + 2*mbc - 2;  
            ifaces_y = my + 2*mbc - 1;
            num_ifaces = ifaces_x*ifaces_y;
        }
        __syncthreads();

        // dtdy0 = mcapa_flag * dtdy/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
                // + (1.0 - mcapa_flag) * dtdy;
        // dtdx0 = mcapa_flag * dtdx/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
                // + (1.0 - mcapa_flag) * dtdx;
        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            int iadd = mbc-1;  // Shift from corner by 1 in each direction
            int I = (iy + iadd)*ys + ix + iadd;

            // int I_capa = I + (mcapa-1)*zs; 
            // // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
            // double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
            // double dtdx_ave = 0.5*(dtdx0 + dtdx_);
            // double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;
            // double dtdy_ave = 0.5*(dtdy0 + dtdy_);

            double dtdy_ave = mcapa_flag * 0.5*(dtdx1d[I-1] + dtdx1d[I]) + (1.0 - mcapa_flag) * dtdx;

            double *const s = start;
            double *const wave = s + mwaves;
            for(int mw = 0; mw < mwaves; mw++)
            {
                /* ------------------------------- Y-directions --------------------------- */
                int I_speeds = I + (mwaves + mw)*zs;
                s[mw] = speeds[I_speeds];

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                    wave[mq] = waves[I_waves];
                }                        

                if (mthlim[mw] > 0)
                {
                    double wnorm2 = 0, dotl = 0, dotr = 0;
                    for(int mq = 0; mq < meqn; mq++)
                    {
                        wnorm2 += pow(wave[mq],2);
                        int I_waves = I + ((mwaves + mw)*meqn + mq)*zs;
                        dotl += wave[mq]*waves[I_waves-ys];
                        dotr += wave[mq]*waves[I_waves+ys];
                    }  
                    wnorm2 = (wnorm2 == 0) ? 1e-15 : wnorm2;  

                    double r = (s[mw] > 0) ? dotl/wnorm2 : dotr/wnorm2;

                    double wlimitr = cudaclaw_limiter(mthlim[mw],r);  

                    for (int mq = 0; mq < meqn; mq++)
                    {
                        wave[mq] *= wlimitr;
                    }                    
                }

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_q = I + mq*zs;
                    double cqyy = (1.0 - fabs(s[mw])*dtdy_ave)*wave[mq];
                    cqyy *= (use_fwaves) ? copysign(1.,s[mw]) : fabs(s[mw]);

                    gm[I_q] += 0.5*cqyy;   
                    gp[I_q] += 0.5*cqyy;  

                    if (order[1] > 0)
                    {                            
                        bmdq_trans[I_q] += cqyy;     
                        bpdq_trans[I_q] -= cqyy;      
                    } 
                }   
            }  
            // dtdy0 = dtdy_; /* update dtdy for next step */
            // dtdx0 = dtdx_; /* update dtdx for next step */
        }  
        __syncthreads();
    }  

    /* ------------------------- First order final update ----------------------------- */

    if (order[1] == 0)
    {

        goto FINAL_UPDATE; /* No transverse propagation; Update the solution and exit */
        // /* No transverse propagation; Update the solution and exit */
        // for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
        // {
        //     int ix = thread_index % mx;
        //     int iy = thread_index/mx;

        //     int iadd = mbc;  // Only update interior cells
        //     int I = (iy + iadd)*ys + (ix + iadd);

        //     int I_capa = I + (mcapa-1)*zs; // mcapa is set to 2 for latlon cordinates (-1 due to the switch between fortran and C)
        //     double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
        //     double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;

        //     for(int mq = 0; mq < meqn; mq++)
        //     {
        //         int I_q = I + mq*zs;
        //         qold[I_q] = qold[I_q] - dtdx_ * (fm[I_q + 1] - fp[I_q])
        //                             - dtdy_ * (gm[I_q + ys] - gp[I_q]);
        //     }        
        // }
        // return;
    }


    /* ------------------------ Transverse Propagation : X-faces ---------------------- */

    // __syncthreads();

    if (threadIdx.x == 0)
    {
        ifaces_x = mx + 2*mbc-1;   /* Visit x - edges of all non-ghost cells */
        ifaces_y = my + 2*mbc-2;                                  
        num_ifaces = ifaces_x*ifaces_y;
    }

    __syncthreads();

    /*     transverse-x

                |     |     | 
                |     |     | 
            ----|-----|-----|-----
                |     X     | 
                |     X  q  |
                |  v--X     |
            ----|--O--|-----|-----
                |     |     |
                |     |     |

    */              
    // dtdx0 = mcapa_flag * dtdx/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    //         + (1.0 - mcapa_flag) * dtdx;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);  /* (ix,iy) = (0,0) maps to first non-ghost value */

        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
        // double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
        // double dtdx_ave = 0.5*(dtdx0 + dtdx_);

        double dtdx_ave = mcapa_flag * 0.5*(dtdx1d[I-1] + dtdx1d[I]) + (1.0 - mcapa_flag) * dtdx;

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const amdq   = ql + meqn;      /* meqn   */
        double *const aux1   = amdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx_ave*bmasdq[mq];
            gm[I_q - 1] -= gupdate;       
            gp[I_q - 1] -= gupdate;   
        }     
        // dtdx0 = dtdx_; /* update dtdx for next step */       
    }

    __syncthreads();

    /*     transverse-x

            |     |     | 
            |     |     | 
        ----|--0--|-----|-----
            |  ^--X     | 
            |     X  q  |
            |     X     |
        ----|-----|-----|-----
            |     |     |
            |     |     |

    */              
    // dtdx0 = mcapa_flag * dtdx/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    //         + (1.0 - mcapa_flag) * dtdx;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);  /* (ix,iy) = (0,0) maps to first non-ghost value */
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
        // double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
        // double dtdx_ave = 0.5*(dtdx0 + dtdx_);

        double dtdx_ave = mcapa_flag * 0.5*(dtdx1d[I-1] + dtdx1d[I]) + (1.0 - mcapa_flag) * dtdx;

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const amdq   = ql + meqn;      /* meqn   */
        double *const aux1   = amdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            amdq[mq] = amdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }         

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx_ave*bpasdq[mq];
            gm[I_q - 1 + ys] -= gupdate;
            gp[I_q - 1 + ys] -= gupdate;
        }     
        // dtdx0 = dtdx_; /* update dtdx for next step */   
    }

    __syncthreads();

    /*     transverse-x

                |     |     | 
                |     |     | 
            ----|-----|-----|-----
                |     X     | 
                |     X  q  |
                |     X--v  |
            ----|-----|--0--|-----
                |     |     |
                |     |     |

        */   
    // dtdx0 = mcapa_flag * dtdx/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    // + (1.0 - mcapa_flag) * dtdx;     

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);  /* (ix,iy) = (0,0) maps to first non-ghost value */
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
        // double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
        // double dtdx_ave = 0.5*(dtdx0 + dtdx_);

        double dtdx_ave = mcapa_flag * 0.5*(dtdx1d[I-1] + dtdx1d[I]) + (1.0 - mcapa_flag) * dtdx;

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const apdq   = ql + meqn;      /* meqn   */
        double *const aux1   = apdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            apdq[mq] = apdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }

            
        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx_ave*bmasdq[mq];
            gm[I_q] -= gupdate;       
            gp[I_q] -= gupdate;
        }
        // dtdx0 = dtdx_; /* update dtdx for next step */
    }

    __syncthreads();

    /*     transverse-x

                |     |     | 
                |     |     | 
            ----|-----|--0--|-----
                |     X--^  | 
                |     X  q  |
                |     X     |
            ----|-----|-----|-----
                |     |     |
                |     |     |

        */              
        // dtdx0 = mcapa_flag * dtdx/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
        // + (1.0 - mcapa_flag) * dtdx;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);  /* (ix,iy) = (0,0) maps to first non-ghost value */
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
        // double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
        // double dtdx_ave = 0.5*(dtdx0 + dtdx_);

        double dtdx_ave = mcapa_flag * 0.5*(dtdx1d[I-1] + dtdx1d[I]) + (1.0 - mcapa_flag) * dtdx;

        double *const qr     = start;          /* meqn   */
        double *const ql     = qr + meqn;      /* meqn   */
        double *const apdq   = ql + meqn;      /* meqn   */
        double *const aux1   = apdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            ql[mq] = qold[I_q-1];
            qr[mq] = qold[I_q];

            apdq[mq] = apdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - ys + (imp - 1)];
                aux2[k] = aux[I_aux      + (imp - 1)];
                aux3[k] = aux[I_aux + ys + (imp - 1)];
            }
        }

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx_ave*bpasdq[mq];
            gm[I_q + ys] -= gupdate;
            gp[I_q + ys] -= gupdate;
        }
        // dtdx0 = dtdx_; /* update dtdx for next step */
        
    } 

    __syncthreads();  


    /* ----------------------------- Transverse : Y-faces ----------------------------- */


    /*  transverse-y

                |     |     
            -----|-----|-----
                |     |     
                |  q  |      
                |     |     
            -----|-XXX-|-----
                | v   |     
                0--   0     
                |     |     
            -----|-----|-----
                |     |     
        */                        

    // __syncthreads();

    if (threadIdx.x == 0)
    {
        ifaces_x = mx + 2*mbc-2;   /* Visit edges of all non-ghost cells */
        ifaces_y = my + 2*mbc-1;
        num_ifaces = ifaces_x*ifaces_y;
    }
    __syncthreads();

    // dtdy0 = mcapa_flag * dtdy/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    //         + (1.0 - mcapa_flag) * dtdy;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;  // Shift from corner by 1 in each direction
        int I =  (iy + iadd)*ys + (ix + iadd);
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;
        // double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;
        // double dtdy_ave = 0.5*(dtdy0 + dtdy_);

        // double dtdy_ave = 0.5*(dtdy1d[I-ys] + dtdy1d[I]);
        double dtdy_ave = mcapa_flag * 0.5*(dtdy1d[I-ys] + dtdy1d[I]) + (1.0 - mcapa_flag) * dtdy;

        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bmdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bmdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy_ave*bmasdq[mq];
            fm[I_q - ys] -= gupdate;        
            fp[I_q - ys] -= gupdate;

        }
        // dtdy0 = dtdy_; /* update dtdy for next step */
    }

    __syncthreads();
    /*  transverse-y

            |     |     
        -----|-----|-----
            |     |     
            |  q  |      
            |     |     
        -----|-XXX-|-----
            |   v |     
            0   --0     
            |     |     
        -----|-----|-----
            |     |     
    */          

    // dtdy0 = mcapa_flag * dtdy/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    //         + (1.0 - mcapa_flag) * dtdy;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;
        // double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;
        // double dtdy_ave = 0.5*(dtdy0 + dtdy_);

        double dtdy_ave = mcapa_flag * 0.5*(dtdy1d[I-ys] + dtdy1d[I]) + (1.0 - mcapa_flag) * dtdy;

        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bmdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bmdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bmdq[mq] = bmdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy_ave*bpasdq[mq];
            fm[I_q - ys + 1] -= gupdate;
            fp[I_q - ys + 1] -= gupdate;                
        }
        // dtdy0 = dtdy_; /* update dtdy for next step */
        
    }

    __syncthreads();

    /*  transverse-y

            |     |     
        -----|-----|-----
            |  q  |     
            O--   |           
            | ^   |     
        -----|-XXX-|-----
            |     |     
            |     | 
            |     |     
        -----|-----|-----
            |     |     
    */        

    // dtdy0 = mcapa_flag * dtdy/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    //         + (1.0 - mcapa_flag) * dtdy;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;
        // double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;
        // double dtdy_ave = 0.5*(dtdy0 + dtdy_);

        double dtdy_ave = mcapa_flag * 0.5*(dtdy1d[I-ys] + dtdy1d[I]) + (1.0 - mcapa_flag) * dtdy;

        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bpdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bpdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bpdq[mq] = bpdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy_ave*bmasdq[mq];
            fm[I_q] -= gupdate;        
            fp[I_q] -= gupdate;
        }
        // dtdy0 = dtdy_; /* update dtdy for next step */
    }

    __syncthreads();

    /*  transverse-y

            |     |     
        -----|-----|-----
            |  q  |     
            |   --0           
            |   ^ |     
        -----|-XXX-|-----
            |     |     
            |     | 
            |     |     
        -----|-----|-----
            |     |     
    */         

    // dtdy0 = mcapa_flag * dtdy/(aux[((0 + (mbc-1))*ys + (0 + (mbc-1))) + (mcapa-1)*zs] + (1.0 - mcapa_flag)) 
    //         + (1.0 - mcapa_flag) * dtdy;
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int iadd = mbc-1;
        int I =  (iy + iadd)*ys + (ix + iadd);
        // int I_capa = I + (mcapa-1)*zs; 
        // // double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;
        // double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;
        // double dtdy_ave = 0.5*(dtdy0 + dtdy_);

        double dtdy_ave = mcapa_flag * 0.5*(dtdy1d[I-ys] + dtdy1d[I]) + (1.0 - mcapa_flag) * dtdy;

        double *const qr     = start;          /* meqn   */
        double *const qd     = qr + meqn;      /* meqn   */
        double *const bpdq   = qd + meqn;      /* meqn   */
        double *const aux1   = bpdq + meqn;    /* 2*maux */
        double *const aux2   = aux1 + 2*maux;  /* 2*maux */
        double *const aux3   = aux2 + 2*maux;  /* 2*maux */
        double *const bmasdq = aux3 + 2*maux;  /* meqn   */
        double *const bpasdq = bmasdq + meqn;  /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;

            qr[mq] = qold[I_q];
            qd[mq] = qold[I_q - ys];

            bpdq[mq] = bpdq_trans[I_q];
        }            

        for(int imp = 0; imp < 2; imp++)
        {
            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                int k = imp*maux + m;
                aux1[k] = aux[I_aux - 1 + ys*(imp - 1)];
                aux2[k] = aux[I_aux     + ys*(imp - 1)];
                aux3[k] = aux[I_aux + 1 + ys*(imp - 1)];
            }
        }

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq, dry_tol);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy_ave*bpasdq[mq];
            fm[I_q + 1] -= gupdate;
            fp[I_q + 1] -= gupdate;
        }   
        // dtdy0 = dtdy_; /* update dtdy for next step */
    } 

    __syncthreads();

    /* ------------------------------- Final update ----------------------------------- */
FINAL_UPDATE: /* No transverse propagation; Update the solution and exit */

    for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
    {
        // Loop over interior cells only
        int ix = thread_index % mx;
        int iy = thread_index/my;

        int iadd = mbc;
        int I = (iy + iadd)*ys + (ix + iadd);

        int I_capa = I + (mcapa-1)*zs; 
        // double dtdx_ = (mcapa > 0) ? dtdx/aux[I_capa] : dtdx;
        // double dtdy_ = (mcapa > 0) ? dtdy/aux[I_capa] : dtdy;
        double dtdx_ = mcapa_flag * dtdx/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdx;
        double dtdy_ = mcapa_flag * dtdy/(aux[I_capa] + (1.0 - mcapa_flag)) + (1.0 - mcapa_flag) * dtdy;

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qold[I_q] = qold[I_q] - dtdx_ * (fm[I_q + 1] - fp[I_q]) 
                                  - dtdy_ * (gm[I_q + ys] - gp[I_q]);

        }        
        //__syncthreads();
#if 1
    if (src2 != NULL && src_term > 0)
    {
        // printf("ix = %d, iy = %d, I = %d\n",ix,iy,I);
        double *const qr = start;          /* meqn   */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qr[mq] = qold[I_q];  
        }
        double *const auxr   = qr + meqn;         /* maux        */
        // for(int m = 0; m < maux; m++)
        // {
        //     /* In case aux is already set */
        //     int I_aux = I + m*zs;
        //     auxr[m] = aux[I_aux];
        // }    

        // First cell in non-ghost cells should be (1,1)
        int i = ix+1;  
        int j = iy+1;
        src2(meqn,maux,xlower,ylower,dx,dy,qr,auxr,t,dt,i,j);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qold[I_q] = qr[mq];  
        }
    }
    // }
#endif
    }

#if 0
    // __syncthreads();

    // if (src2 != NULL && src_term > 0)
    // {
    //     for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
    //     {
    //         // Loop over interior cells only
    //         int ix = thread_index % mx;
    //         int iy = thread_index/my;

    //         int iadd = mbc;
    //         int I = (iy + iadd)*ys + (ix + iadd);
            
    //             // printf("ix = %d, iy = %d, I = %d\n",ix,iy,I);
    //             double *const qr = start;          /* meqn   */
    //             for(int mq = 0; mq < meqn; mq++)
    //             {
    //                 int I_q = I + mq*zs;
    //                 qr[mq] = qold[I_q];  
    //             }
    //             double *const auxr   = qr + meqn;         /* maux        */
    //             for(int m = 0; m < maux; m++)
    //             {
    //                 /* In case aux is already set */
    //                 int I_aux = I + m*zs;
    //                 auxr[m] = aux[I_aux];
    //             }    

    //             // First cell in non-ghost cells should be (1,1)
    //             int i = ix+1;  
    //             int j = iy+1;
    //             src2(meqn,maux,xlower,ylower,dx,dy,qr,auxr,t,dt,i,j);

    //             for(int mq = 0; mq < meqn; mq++)
    //             {
    //                 int I_q = I + mq*zs;
    //                 qold[I_q] = qr[mq];  
    //             }
    //         }
    // }
#endif
}

/* ---------------------------------------------------------------------------------------
   PUBLIC function  
   ------------------------------------------------------------------------------------ */
__global__
void cudaclaw_flux2_and_update_batch (const int mx,    const int my, 
                                      const int meqn,  const int mbc, 
                                      const int maux,  const int mwaves, 
                                      const int mwork, const double dt, 
                                      const double t,  const int src_term,
                                      const int mcapa, const double dry_tol,
                                      cudaclaw_fluxes_t* array_fluxes_struct,
                                      double * maxcflblocks,
                                      cudaclaw_cuda_rpn2_t rpn2,
                                      cudaclaw_cuda_rpt2_t rpt2,
                                      cudaclaw_cuda_b4step2_t b4step2,
                                      cudaclaw_cuda_src2_t src2)
    {
        cudaclaw_flux2_and_update(mx,my,meqn,mbc,maux,mwaves,mwork,
                                  src_term, mcapa, dry_tol,
                                  array_fluxes_struct[blockIdx.z].xlower,
                                  array_fluxes_struct[blockIdx.z].ylower,
                                  array_fluxes_struct[blockIdx.z].dx,
                                  array_fluxes_struct[blockIdx.z].dy,
                                  array_fluxes_struct[blockIdx.z].qold_dev,
                                  array_fluxes_struct[blockIdx.z].aux_dev,
                                  array_fluxes_struct[blockIdx.z].fm_dev,
                                  array_fluxes_struct[blockIdx.z].fp_dev,
                                  array_fluxes_struct[blockIdx.z].gm_dev,
                                  array_fluxes_struct[blockIdx.z].gp_dev,
                                  array_fluxes_struct[blockIdx.z].amdq_dev,
                                  array_fluxes_struct[blockIdx.z].apdq_dev,
                                  array_fluxes_struct[blockIdx.z].bmdq_dev,
                                  array_fluxes_struct[blockIdx.z].bpdq_dev,
                                  array_fluxes_struct[blockIdx.z].waves_dev,
                                  array_fluxes_struct[blockIdx.z].speeds_dev, 
                                  array_fluxes_struct[blockIdx.z].dtdx1d_dev,
                                  array_fluxes_struct[blockIdx.z].dtdy1d_dev,
                                  maxcflblocks,
                                  rpn2, rpt2, b4step2, src2, t, dt);
}

