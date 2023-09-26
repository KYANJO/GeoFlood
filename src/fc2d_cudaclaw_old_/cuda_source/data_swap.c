/* 
@author: Brian Kyanjo
@date: 2023/15/08
@purpose: This program reads in a (m,i,j) or (i,j,m) array and swaps the indices
          to the other form. This is useful for exchanging data between geoclaw (m,i,j) 
          routines and cudaclaw (i,j,m)routines
*/

void swap_mij2ijm(int mx, int my, int mbc, int meqn, int maux, double qold[meqn][mx+2*mbc][my+2*mbc],
                double qold_transpose[mx+2*mbc][my+2*mbc][meqn], double aux[maux][mx+2*mbc][my+2*mbc],
               double aux_transpose[mx+2*mbc][my+2*mbc][maux])
{
    int i,j,m;

    // swap qold
        
    for(i=0; i<mx+2*mbc;i++)
    {
        for(j=0; j<my+2*mbc;j++)
        {   
            for(m=0;m<meqn;m++)
            {
                qold_transpose[i][j][m] = qold[m][i][j];
            }
        }
    }

    // swap aux
    for(i=0; i<mx+2*mbc;i++)
    {
        for(j=0; j<my+2*mbc;j++)
        {   
            for(m=0;m<maux;m++)
            {
                aux_transpose[i][j][m] = aux[m][i][j];
            }
        }
    }

}

void swap_ijm2mij(int mx, int my, int mbc, int meqn, int maux, double qold[mx+2*mbc][my+2*mbc][meqn],
                double qold_transpose[meqn][mx+2*mbc][my+2*mbc], double aux[mx+2*mbc][my+2*mbc][maux],
                double aux_transpose[maux][mx+2*mbc][my+2*mbc])
{
    int i,j,m;

    // swap qold
    for(m=0;m<meqn;m++)
    {
        for(i=0; i<mx+2*mbc;i++)
        {
            for(j=0; j<my+2*mbc;j++)
            {    
                qold_transpose[m][i][j] = qold[i][j][m];
            }
        }
    }
    
    // swap aux
    for(m=0;m<maux;m++)
    {
        for(i=0; i<mx+2*mbc;i++)
        {
            for(j=0; j<my+2*mbc;j++)
            {  
                aux_transpose[m][i][j] = aux[i][j][m];
            }
        }
    }

}
