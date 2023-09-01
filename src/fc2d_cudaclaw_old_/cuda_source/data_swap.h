// call swap functions
#ifndef DATA_SWAP_H
#define DATA_SWAP_H

#ifdef __cplusplus
extern "C" 
{
#if 0
}
#endif
#endif

void swap_ijm2mij(int mx, int my, int mbc, int meqn, int maux, double qold[mx+2*mbc][my+2*mbc][meqn],
                double qold_transpose[meqn][mx+2*mbc][my+2*mbc], double aux[mx+2*mbc][my+2*mbc][maux],
                double aux_transpose[maux][mx+2*mbc][my+2*mbc]);

void swap_mij2ijm(int mx, int my, int mbc, int meqn, int maux, double qold[meqn][mx+2*mbc][my+2*mbc],
                double qold_transpose[mx+2*mbc][my+2*mbc][meqn], double aux[maux][mx+2*mbc][my+2*mbc],
               double aux_transpose[mx+2*mbc][my+2*mbc][maux]);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif
#endif // DATA_SWAP_H