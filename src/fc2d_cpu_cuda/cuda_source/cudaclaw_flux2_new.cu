#include "../fc2d_cudaclaw_cuda.h"

#include "../fc2d_cudaclaw_check.h"

#include <cub/block/block_reduce.cuh>  

#include "cudaclaw_allocate.h"  /* Needed to for definition of 'fluxes' */

__constant__ int order[2];
__constant__ int mthlim[FC2D_CUDACLAW_MWAVES];
__constant__ int use_fwaves;

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
                               const int mwork,
                               const double xlower, const double ylower, 
                               const double dx,     const double dy,
                               double *const qold,       double *const aux, 
                               double *const fm,         double *const fp, 
                               double *const gm,         double *const gp,
                               double *const amdq_trans, double *const apdq_trans, 
                               double *const bmdq_trans, double *const bpdq_trans,
                               double *const waves,      double *const speeds,
                               double *const maxcflblocks,
                               cudaclaw_cuda_rpn2_t rpn2,
                               cudaclaw_cuda_rpt2_t rpt2,
                               cudaclaw_cuda_b4step2_t b4step2,
                               double t,double dt)
{
    typedef cub::BlockReduce<double,FC2D_CUDACLAW_BLOCK_SIZE> BlockReduce;
    
    __shared__ typename BlockReduce::TempStorage temp_storage;

    extern __shared__ double shared_mem[];

    double* start  = shared_mem + mwork*threadIdx.x;

    /* --------------------------------- Start code ----------------------------------- */

    __shared__ double dtdx, dtdy;
    __shared__ int xs,ys,zs;
    __shared__ int ifaces_x, ifaces_y, num_ifaces;

    if (threadIdx.x == 0)
    {
        dtdx = dt/dx;
        dtdy = dt/dy;        

        /* Compute strides */
        xs = 1;
        ys = (2*mbc + mx)*xs;
        zs = (2*mbc + my)*xs*ys;

        ifaces_x = mx + 2*mbc-1;
        ifaces_y = my + 2*mbc-1;
        num_ifaces = ifaces_x*ifaces_y;

    }
   
    __syncthreads();


    double maxcfl = 0;
    double *const qr     = start; 
    double *const auxr   = qr      + meqn;         /* maux        */
    double *const ql     = auxr   + maux;         /* meqn        */
    double *const auxl   = ql     + meqn;         /* maux        */
    double *const s      = auxl   + maux;         /* mwaves      */
    double *const wave   = s      + mwaves;       /* meqn*mwaves */
    double *const amdq   = wave   + meqn*mwaves;  /* meqn        */
    double *const apdq   = amdq   + meqn;    

    /* -------------------------- Compute fluctuations -------------------------------- */
    // for(int mq = 0; mq < meqn; mq++)
    // {
    //     amdq[mq] = 0.0;
    //     apdq[mq] = 0.0;
    // }

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    {
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I = (iy + 1)*ys + (ix + 1);  /* Start one cell from left/bottom edge */

        // // // double *const qr     = start;                 /* meqn        */
        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qr[mq] = qold[I_q];        /* Right */
        }

        // double *const auxr   = qr      + meqn;         /* maux        */
        for(int m = 0; m < maux; m++)
        {
            int I_aux = I + m*zs;
            auxr[m] = aux[I_aux];
        }               

        {
            /* ------------------------ Normal solve in X direction ------------------- */
            // double *const ql     = auxr   + maux;         /* meqn        */
            // double *const auxl   = ql     + meqn;         /* maux        */
            // double *const s      = auxl   + maux;         /* mwaves      */
            // double *const wave   = s      + mwaves;       /* meqn*mwaves */
            // double *const amdq   = wave   + meqn*mwaves;  /* meqn        */
            // double *const apdq   = amdq   + meqn;         /* meqn        */
            // int I = (iy + 1 )*ys + (ix + 1); 
            
            // if (ix < mx + 2 && iy < my + 1){

                for(int mq = 0; mq < meqn; mq++)
                {
                    int I_q = I + mq*zs;
                    ql[mq] = qold[I_q - 1];    /* Left  */
                }

                for(int m = 0; m < maux; m++)
                {
                    int I_aux = I + m*zs;
                    auxl[m] = aux[I_aux-1];
                }               
            
                rpn2(0, meqn, mwaves, maux, ql, qr, auxl, auxr, wave, s, amdq, apdq,ix,iy);

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
            

                for(int mw = 0; mw < mwaves; mw++)
                {
                    if (ix > 0 && ix < mx + 1)
                    {
                        maxcfl = max(maxcfl,fabs(s[mw]*dtdx));
                    }
                    // maxcfl = max(maxcfl,fabs(s[mw]*dtdx));

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
               
            // }
        }

        {

            /* ------------------------ Normal solve in Y direction ------------------- */
            double *const qd     = auxr   + maux;         /* meqn        */
            double *const auxd   = qd     + meqn;         /* maux        */
            double *const s      = auxd   + maux;         /* mwaves      */
            double *const wave   = s      + mwaves;       /* meqn*mwaves */
            double *const bmdq   = wave   + meqn*mwaves;  /* meqn        */
            double *const bpdq   = bmdq   + meqn;         /* meqn        */

            int I = (iy + 1)*ys + (ix + 0)*xs; 

            for(int mq = 0; mq < meqn; mq++)
            {
                int I_q = I + mq*zs;
                qd[mq] = qold[I_q - ys];   /* Down  */  
            }

            for(int m = 0; m < maux; m++)
            {
                int I_aux = I + m*zs;
                auxd[m] = aux[I_aux - ys];
            }      
            
            rpn2(1, meqn, mwaves, maux, qd, qr, auxd, auxr, wave, s, bmdq, bpdq,ix,iy);

            /* Set value at bottom interface of cell I */
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

            for(int mw = 0; mw < mwaves; mw++)
            {
                if (ix > 0 && ix < mx + 1)
                {
                    maxcfl = max(maxcfl,fabs(s[mw]*dtdy));
                }
                // maxcfl = max(maxcfl,fabs(s[mw])*dtdy);

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
        }
    }
 
    maxcflblocks[blockIdx.z] = BlockReduce(temp_storage).Reduce(maxcfl,cub::Max());
    // printf("maxcflblocks[blockIdx.z] = %f, maxcfl = %f, dt = %f \n",maxcflblocks[blockIdx.z],maxcfl,dt);
    /* ---------------------- Second order corrections and limiters --------------------*/  

    if (order[0] == 2)
    {
        if (threadIdx.x == 0){
            ifaces_x = mx + 1;  
            ifaces_y = my + 2;
            num_ifaces = ifaces_x*ifaces_y;
        }
        __syncthreads();

        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            int I = (iy + mbc-1)*ys + (ix + mbc);

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
                    double cqxx = (1.0 - fabs(s[mw])*dtdx)*wave[mq];
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
        }
    }

    __syncthreads();

    if (order[0] == 2)
    {
        if (threadIdx.x == 0){
            ifaces_x = mx + 2;  
            ifaces_y = my + 1;
            num_ifaces = ifaces_x*ifaces_y;
        }
        __syncthreads();

        for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
        { 
            int ix = thread_index % ifaces_x;
            int iy = thread_index/ifaces_x;

            /* Start at first non-ghost interior cell */
            int I = (iy + mbc)*ys + ix + mbc - 1;

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
                    double cqyy = (1.0 - fabs(s[mw])*dtdx)*wave[mq];
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
        }  
        __syncthreads();
    }  


    /* ------------------------- First order final update ----------------------------- */

    if (order[1] == 0)
    {
        /* No transverse propagation; Update the solution and exit */
        for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
        {
            int ix = thread_index % mx;
            int iy = thread_index/mx;

            int I = (ix + mbc) + (iy + mbc)*ys;

            for(int mq = 0; mq < meqn; mq++)
            {
                int I_q = I + mq*zs;
                qold[I_q] = qold[I_q] - dtdx * (fm[I_q + 1] - fp[I_q]) 
                                      - dtdy * (gm[I_q + ys] - gp[I_q]); 
            }        
        }
        return;
    }



    /* ------------------------ Transverse Propagation : X-faces ---------------------- */
    if (threadIdx.x == 0){
        ifaces_x = mx + 1;   /* Visit x - edges of all non-ghost cells */
        ifaces_y = my + 2;
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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

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

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bmasdq[mq];

            gm[I_q - 1] -= gupdate;       
            gp[I_q - 1] -= gupdate;   

            /* up and down update (left hand side) */
            // atomicAdd(&gm[I_q-1],-gupdate);
            // atomicAdd(&gp[I_q-1],-gupdate);
 
            // gupdate = 0.5*dtdx*bpasdq[mq];
            // atomicAdd(&gm[I_q - 1 + ys],-gupdate);
            // atomicAdd(&gp[I_q - 1 + ys],-gupdate);
        }            
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
  /* this for loop is not needed if we use atomicAdd for the up and down update */
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

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

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,0,amdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bpasdq[mq];
            gm[I_q - 1 + ys] -= gupdate;
            gp[I_q - 1 + ys] -= gupdate;
        }        
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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

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

             
        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bmasdq[mq];
            gm[I_q] -= gupdate;       
            gp[I_q] -= gupdate;

            /* up and down update (right hand side) */
            // atomicAdd(&gm[I_q],-gupdate);       
            // atomicAdd(&gp[I_q],-gupdate);

            // gupdate = 0.5*dtdx*bpasdq[mq];
            // atomicAdd(&gm[I_q + ys],-gupdate);
            // atomicAdd(&gp[I_q + ys],-gupdate);

        }
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
/*--  this for loop is not needed if we use atomicAdd for the up and down update -- */
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc-1)*ys + (ix + mbc);  /* (ix,iy) = (0,0) maps to first non-ghost value */

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

        rpt2(0,meqn,mwaves,maux,ql,qr,aux1,aux2,aux3,1,apdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdx*bpasdq[mq];
            gm[I_q + ys] -= gupdate;
            gp[I_q + ys] -= gupdate;
        }
        
    } 

    /* May not the synchthreads(), below, since gm/gp updated above, but only fm/gp
       updated below */
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

    if (threadIdx.x == 0){
        ifaces_x = mx + 2;   /* Visit edges of all non-ghost cells */
        ifaces_y = my + 1;
        num_ifaces = ifaces_x*ifaces_y;
    }
    __syncthreads();

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);

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

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bmasdq[mq];
            fm[I_q - ys] -= gupdate;        
            fp[I_q - ys] -= gupdate;

            /* left and right update down */
            // atomicAdd(&fm[I_q - ys],-gupdate);        
            // atomicAdd(&fp[I_q - ys],-gupdate);

            // gupdate = 0.5*dtdy*bpasdq[mq];
            // atomicAdd(&fm[I_q - ys + 1],-gupdate);
            // atomicAdd(&fp[I_q - ys + 1],-gupdate);   
        }

    }

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

    /* This for loop is not needed if we use atomicAdd for the left and right update down */
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


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

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,0,bmdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bpasdq[mq];
            fm[I_q - ys + 1] -= gupdate;
            fp[I_q - ys + 1] -= gupdate;                
        }
        
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

    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


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

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bmasdq[mq];
            fm[I_q] -= gupdate;        
            fp[I_q] -= gupdate;

            /* left and right update up */
            // // atomicAdd(&fm[I_q], -gupdate);        
            // // atomicAdd(&fp[I_q], -gupdate);

            // gupdate = 0.5*dtdy*bpasdq[mq];
            // // atomicAdd(&fm[I_q + 1],-gupdate);
            // // atomicAdd(&fp[I_q + 1],-gupdate);
            // fm[I_q+1] -= gupdate;        
            // fp[I_q+1] -= gupdate;
        }
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

    /* This for loop is not needed if we use atomicAdd for the left and right update up */
    for(int thread_index = threadIdx.x; thread_index < num_ifaces; thread_index += blockDim.x)
    { 
        int ix = thread_index % ifaces_x;
        int iy = thread_index/ifaces_x;

        int I =  (iy + mbc)*ys + (ix + mbc-1);


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

        rpt2(1,meqn,mwaves,maux,qd,qr,aux1,aux2,aux3,1,bpdq,bmasdq,bpasdq,ix,iy);

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;  
            double gupdate = 0.5*dtdy*bpasdq[mq];
            fm[I_q + 1] -= gupdate;
            fp[I_q + 1] -= gupdate;
        }   
        
    } 

    __syncthreads();

    /* ------------------------------- Final update ----------------------------------- */

    for(int thread_index = threadIdx.x; thread_index < mx*my; thread_index += blockDim.x)
    {
        int ix = thread_index % mx;
        int iy = thread_index/my;

        int I = (ix + mbc)*xs + (iy + mbc)*ys;

        for(int mq = 0; mq < meqn; mq++)
        {
            int I_q = I + mq*zs;
            qold[I_q] = qold[I_q] - dtdx * (fm[I_q + 1] - fp[I_q]) 
                                  - dtdy * (gm[I_q + ys] - gp[I_q]);
        }        
    }

}


/* ---------------------------------------------------------------------------------------
   PUBLIC function  
   ------------------------------------------------------------------------------------ */
__global__
void cudaclaw_flux2_and_update_batch (const int mx,    const int my, 
                                        const int meqn,  const int mbc, 
                                        const int maux,  const int mwaves, 
                                        const int mwork,
                                        const double dt, const double t,
                                        cudaclaw_fluxes_t* array_fluxes_struct,
                                        double * maxcflblocks,
                                        cudaclaw_cuda_rpn2_t rpn2,
                                        cudaclaw_cuda_rpt2_t rpt2,
                                        cudaclaw_cuda_b4step2_t b4step2)
    {
        cudaclaw_flux2_and_update(mx,my,meqn,mbc,maux,mwaves,mwork,
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
                                    maxcflblocks,
                                    rpn2, rpt2, b4step2, t,dt);
}



