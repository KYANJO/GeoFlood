/*
@author: Brian Kyanjo
@description: Solves a shallow water equation (swe) given a single left and right states
@date: 30th July 2023
@reference: Solver is described in J. Comput.Phys. (6): 3089-3113, March 2008 Augmented 
            Riemann Solvers for the swe with steady states and Inundation

@note: - To use the original solver call with maxiter=1.
       - This solver allows iteration when maxiter > 1. The iteration seems to help with instabilities that arise (with any solver) as flow becomes transcritical over variable topography due to loss of hyperbolicity.
*/

// void riemann_aug_JCP(int maxiter, int meqn, )










// =============================================================================
void riemanntype(double hL, double hR, double uL, double uR, double hm, 
                 double s1m, double s2m, bool rare1, bool rare2, int maxiter,
                 double drytol, double g)
{
    // Local variables
    double u1m,u2m,delu;
    double h_max, h_min, h0, F_max, F_min, dfdh, F0, slope, gL, gR;
    int iter;

    // Test for Riemann structure
    h_min = fmin(hR, hL);
    h_max = fmax(hR, hL);
    delu = uR - uL;

    if (h_min <= drytol) 
    {
        hm = 0.0;
        um = 0.0;
        s1m = uR + uL - 2.0*sqrt(g*hR) + 2.0*sqrt(g*hL);
        s2m = uR + uL - 2.0*sqrt(g*hR) - 2.0*sqrt(g*hL);
        if (hL<=0.0)
        {
            rare2 = true;
            rare1 = false;
        }
        else
        {
            rare1 = true;
            rare2 = false;
        }
    }
    else
    {
        F_min = delu + 2.0*(sqrt(g*h_min) - sqrt(g*h_max));
        F_max = delu + (h_max - h_min)*(sqrt(0.5*g*(h_max + h_min)/(h_max*h_min)));

        if (F_min > 0.0) //2-rarefactions
        {
            hm = (1/(16*g))*pow(fmax(0.0,-delu+2*(sqrt(g*hL)+srt(g*hR))),2);
            double sign_hm = (hm >= 0.0) ? 1.0 : -1.0;
            um = sign_hm*(uL+2*(sqrt(g*hL)-sqrt(g*hm)));

            s1m = uL + 2*sqrt(g*hL) - 3*sqrt(g*hm);
            s2m = uR - 2*sqrt(g*hR) + 3*sqrt(g*hm);

            rare1 = true;
            rare2 = true;
        }
        else if (F_max <= 0.0) // 2-shocks
        {
            // root finding using a Newton iteration on sqrt(h)
            h0 = h_max;
            for (iter=0; iter < maxiter; ++iter)
            {
                gL = sqrt(0.5*g*(1/h0 + 1/hL));
                gR = sqrt(0.5*g*(1/h0 + 1/hR));
                F0 = delu + (h0 - hL)*gL + (h0 - hR)*gR;
                dfdh = gL - g*(h0 - hL)/(4.0*h0*h0*gL) + gR - g*(h0 - hR)/(4.0*h0*h0*gR);
                slope = 2.0*sqrt(h0)*dfdh;
                h0 = pow((sqrt(h0) - F0/slope),2);
            }
            hm = h0;
            u1m = uL - (hm-hL)*sqrt((0.5*g)*(1/hm + 1/hL));
            u2m = uR + (hm-hR)*sqrt((0.5*g)*(1/hm + 1/hR));
            um = 0.5*(u1m + u2m);

            s1m = u1m - 2*sqrt(g*hm);
            s2m = u2m + 2*sqrt(g*hm);

            rare1 = false;
            rare2 = false;
        }
        else // 1-shock 1-rarefaction
        {
            h0 = h_min;
            for(iter=0; iter < maxiter; ++iter)
            {
                F0 = delu + 2.0*(sqrt(g*h0) - sqrt(g*h_max)) + (h0-h_min)*sqrt(0.5*g*(1/h0 + 1/h_min));
                slope = (F_max - F0)/(h_max - h_min);
                h0 = h0 - F0/slope;
            }
            hm = h0;
            if (hL > hR)
            {
                um = uL + 2.0*(sqrt(g*hL) - sqrt(g*hm));
                s1m = uL + 2.0*sqrt(g*hL) - 3.0*sqrt(g*hm);
                s2m = uL + 2.0*sqrt(g*hL) - sqrt(g*hm);
                rare1 = true;
                rare2 = false;
            }
            else
            {
                s2m = uR - 2.0*sqrt(g*hR) + 3.0*sqrt(g*hm);
                s1m = uR - 2.0*sqrt(g*hR) + sqrt(g*hm);
                um = uR - 2.0*(sqrt(g*hR) + 2.0*sqrt(g*hm));
                rare1 = false;
                rare2 = true;
            }
        }
    }
}