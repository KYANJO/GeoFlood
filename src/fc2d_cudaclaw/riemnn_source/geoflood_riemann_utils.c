/*
@author: Brian Kyanjo
@description: Solves a shallow water equation (swe) given a single left and right states
@date: 30th July 2023
@reference: Solver is described in J. Comput.Phys. (6): 3089-3113, March 2008 Augmented 
            Riemann Solvers for the swe with steady states and Inundation
*/

/* === Begin fuction riemann_aug_JCP======================================================== @description: - Solves swe give single left and right states
@note: - To use the original solver call with maxiter=1.
       - This solver allows iteration when maxiter > 1. The iteration seems to help  
         with instabilities that arise (with any solver) as flow becomes transcritical 
         over variable topography due to loss of hyperbolicity. 
*/
void riemann_aug_JCP(int maxiter, int meqn, int mwaves, double hL, double hR, double huL, 
                     double huR, double hvL, double hvR, double bL, double bR, double uL, 
                     double uR, double vL, double vR, double phiL, double phiR, double sE1, double sE2, double drytol, double g, double sw[], double fw[][])
{
    // Local variables
    integer m, mw, k, iter;
    double A[2][2];
    double r[2][2];
    double lambda[2];
    double del[2];
    double beta[2];

    double delh,delhu,delphi,delb,delnorm;
    double rare1st,rare2st,sdelta,raremin,raremax;
    double criticaltol,convergencetol,raretol;
    double s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar;
    double huRstar,huLstar,uRstar,uLstar,hstarHLL;
    double deldelh,deldelphi;
    double s1m,s2m,hm;
    double det1,det2,det3,determinant;

    bool rare1,rare2,rarecorrector,rarecorrectortest,sonic;

    // determine del vectors
    delh = hR - hL;
    delhu = huR - huL;
    delphi = phiR - phiL;
    delb = bR - bL;
    delnorm = delh*delh + delphi*delphi;

    //  call riemanntype function
    riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,0,drytol,g)

    lambda[0] = fmin(sE1,s2m); // Modified Einfeldt speed
    lambda[2] = fmax(sE2,s1m); // Modified Einfeldt speed
    sE1 = lambda[0];
    sE2 = lambda[2];
    lambda[1] = 0.0 // Fix to avoid uninitialized value in loop on mw

    hstarHLL = fmax((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.0); // middle state in a HLL solver

    // determine the middle entropy corrector wave
    rarecorrectortest = false;
    rarecorrector = false;
    if (rarecorrectortest)
    {
        sdelta = lambda[2] - lambda[0];
        raremin = 0.5;
        raremax = 0.9;
        if (rare1 && sE1*s1m < 0.0) raremin = 0.2;
        if (rare2 && sE2*s2m < 0.0) raremin = 0.2;
        if (rare1 || rare2)
        {
            // see which rarefaction is larger
            rare1st = 3.0*(sqrt(g*hL) - sqrt(g*hm));
            rare2st = 3.0*(sqrt(g*hR) - sqrt(g*hm));
            if (fmax(rare1st,rare2st) > raremin*sdelta && fmax(rare1st,rare2st) < raremax*sdelta)
            {
                rarecorrector = true;
                if (rare1st > rare2st)
                {
                    lambda[1] = s1m;
                }
                else if (rare2st > rare1st)
                {
                    lambda[1] = s2m;
                }
                else
                {
                    lambda[1] = 0.5*(s1m + s2m);
                }
            }
        }
        if (hstarHLL < fmin(hL,hR)/5.0) rarecorrector = false;
    }

    for (mw = 0; mw<=mwaves; ++mw)
    {
        r[0][mw] = 1.0;
        r[1][mw] = lambda[mw];
        r[2][mw] = pow(lambda[mw],2);
    }

    if (!rarecorrector)
    {
        lambda[1]= 0.5*(lambda[0] + lambda[2]);
        r[0][1] = 0.0;
        r[1][1] = 0.0;
        r[2][1] = 1.0;
    }

    // determine the steady state wave
    criticaltol = 1.0e-6;
    deldelh = -delh;
    deldelphi = -g*0.5*(hR+hL)*delb;

    // determine a few quantities needed for steady state wave if iterated
    hLstar = hL;
    hRstar = hR;
    uLstar = uL;
    uRstar = uR;
    huLstar = uLstar*hLstar;
    huRstar = uRstar*hRstar;

    // iterate to better determine steady state wave
    convergencetol = 1.0e-6;
    for (iter=0; iter<=maiter; ++iter)
    {
        // determine steady state wave (this will be subtracted from the delta vectors)
        if (fmin(hLstar,hRstar) < drytol && rarecorrector)
        {
            rarecorrector = false;
            hLstar = hL;
            hRstar = hR;
            uLstar = uL;
            uRstar = uR;
            huLstar = uLstar*hLstar;
            huRstar = uRstar*hRstar;
            lambda[1] = 0.5*(lambda[0] + lambda[2]);
            r[0][1] = 0.0;
            r[1][1] = 0.0;
            r[2][1] = 1.0;
        }

        hbar = fmax(0.5*(hLstar+hRstar),0.0);
        s1s2bar = 0.25*pow((uLstar+uRstar),2) - g*hbar;
        s1s2tilde = fmax(0.0,uLstar*uRstar) - g*hbar;

        //  find if sonic problem
        sonic = false;
        if (fabs(s1s2bar) <= criticaltol) sonic = true;
        if (s1s2bar*s1s2tilde <= criticaltol) sonic = true;
        if (s1s2bar*sE1*sE2 <= criticaltol) sonic = true;
        if (fmin(fabs(sE1),fabs(sE2)) <= criticaltol) sonic = true;
        if (sE1 < 0.0 && s1m > 0.0) sonic = true;
        if (sE2 > 0.0 && s2m < 0.0) sonic = true;
        if ((uL+sqrt(g*hL))*(uR+sqrt(g*hR)) < 0.0) sonic = true;
        if ((uL-sqrt(g*hL))*(uR-sqrt(g*hR)) < 0.0) sonic = true;

        // find jump in h, deldelh
        if (sonic)
        {
            deldelh = -delh;
        }
        else
        {
            deldelh = delb*g*hbar/s1s2bar;
        }

        //  find bounds in case of critical state resonance, or negative states
        if (sE1 < -criticaltol && sE2 > criticaltol)
        {
            deldelh = fmin(deldelh,hstarHLL*(sE2 - sE1)/sE2);
            deldelh = fmax(deldelh,hstarHLL*(sE2 - sE1)/sE1);
        }
        else if (sE1 >= criticaltol)
        {
            deldelh = fmin(deldelh,hstarHLL*(sE2 - sE1)/sE1);
            deldelh = fmax(deldelh,-hL);
        }
        else if (sE2 <= -criticaltol)
        {
            deldelh = fmin(deldelh,hR);
            deldelh = fmax(deldelh,hstarHLL*(sE2 - sE1)/sE2);
        }

        //  find jump in phi, deldelphi
        if (sonic)
        {
            deldelphi = -g*hbar*delb;
        }
        else
        {
            deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar;
        }

        // find bounds in case of critical state resonance, or negative states
        deldelphi = fmin(deldelphi,g*fmax(-hLstar*delb,-hRstar*delb));
        deldelphi = fmax(deldelphi,g*fmin(-hLstar*delb,-hRstar*delb));

        del[0] = delh - deldelh;
        del[1] = delhu;
        del[2] = delphi - deldelphi;

        //--- determine determinant of eigen matrix ----
        det1 = r[0][0]*(r[1][1]*r[2][2] - r[1][2]*r[2][1]);
        det2 = r[0][1]*(r[1][0]*r[2][2] - r[1][2]*r[2][0]);
        det3 = r[0][2]*(r[1][0]*r[2][1] - r[1][1]*r[2][0]);
        determinant = det1 - det2 + det3;

        // --- solve for beta(k) using Cramers Rule ---
        for(k=0; k<=2; ++k)
        {
            for(mw=0; mw<=2; ++mw)
            {
                A[0][mw] = r[0][mw];
                A[1][mw] = r[1][mw];
                A[2][mw] = r[2][mw];
            }
            A[0][k] = del[0];
            A[1][k] = del[1];
            A[2][k] = del[2];
            det1 = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]);
            det2 = A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
            det3 = A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
            beta[k] = (det1 - det2 + det3)/determinant;
        }

        //  exit if things aren't changing
        if (fabs(pow(del[0],2)+pow(del[2],2) - delnorm) < convergencetol) break;

        delnorm = pow(del[0],2)+pow(del[2],2);

        // find new states qLstar and qRstar on either side of the interface
        hLstar = hL;
        hRstar = hR;
        uLstar = uL;
        uRstar = uR;
        huLstar = uLstar*hLstar;
        huRstar = uRstar*hRstar;
        for (mw=0; mw<=mwaves; ++mw)
        {
            if (lambda[mw] < 0.0)
            {
                hLstar = hLstar + beta[mw]*r[0][mw];
                huLstar = huLstar + beta[mw]*r[1][mw];
            }
        }

        for (mw = mwaves-1; mw >= 0; mw--)
        {
            if (lambda[mq] > 0.0)
            {
                hRstar = hRstar - beta[mw]*r[0][mw];
                huRstar = huRstar - beta[mw]*r[1][mw];
            }
        }

        if (hLstar > drytol) 
        {
            uLstar = huLstar/hLstar;
        }
        else
        {
            hLstar = fmax(hLstar,0.0);
            uLstar = 0.0;
        }

        if (hRstar > drytol) 
        {
            uRstar = huRstar/hRstar;
        }
        else
        {
            hRstar = fmax(hRstar,0.0);
            uRstar = 0.0;
        }
    } // --end iteration on Riemann problem

    for (mw=0; mw<=mwaves; ++mw)
    {
        sw[mw] = lambda[mw];
        fw[0,mw] = beta[mw]*r[1][mw];
        fw[1,mw] = beta[mw]*r[2][mw];
        fw[2,mw] = beta[mw]*r[1][mw];
    }

    // find transverse components (ie huv jumps)
    fw[2][0] = fw[2][0]*vL;
    fw[2][2] = fw[2][2]*vR;
    fw[2][2] = hR*uR*vR - hL*uL*vL - fw[2][0] - fw[2][2];
}// === end function riemann_aug_JCP ========================================================


// === Begin fuction Riemann_ssqfwave========================================================
// @description: - Solves swe give single left and right states
//               - Steady state wave is subtracted from delta [q,f]^T before decompositon
void riemann_ssqfwave(int maxiter, int meqn, int mwaves, double hL, double hR, double huL, 
                      double huR, double hvL, double hvR, double bL, double bR, double uL, 
                      double uR, double vL, double vR, double phiL, double phiR, double sE1, double sE2, double drytol, double g, double sw[], double fw[][])
{
    // Local variables
    double delh, delhu, delphi, delb, delhdecomp, delphidecomp, deldelh;
    double s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar;
    double uRstar,uLstar,hstarHLL;
    double deldelh, deldelphi;
    double alpha1, alpha2, beta1, beta2,delalpha1,delalpha2;
    double criticaltol, convergencetol;
    double sL,sR;
    double uhat,chat,sRoe1,sRoe2;

    integer iter;
    bool sonic;

    // determine del vectors
    delh = hR - hL;
    delhu = huR - huL;
    delphi = phiR - phiL;
    delb = bR - bL;

    convergencetol = 1.0e-16;
    criticaltol = 1.0e-99;

    deldelh = -delb;
    deldelphi = -g*0.5*(hR + hL)*delb;

    // if no source term, skip determining steady state wave
    if (fabs(delb) > 0.0)
    {
        // determine a few quantities needed for steady state wave if iterated
        hLstar = hL;
        hRstar = hR;
        uLstar = uL;
        uRstar = uR;
        hstarHLL = fmax((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.0); // middle state in a HLL solver

        alpha1 = 0.0;
        alpha2 = 0.0;

        // iterate to better determine the Riemann problem
        for (iter=0; iter <= maxiter; ++iter)
        {
            // determine steady state wave (this will be subtracted from delta vectors)
            hbar = fmax(0.5*(hLstar + hRstar),0.0);
            s1s2bar = 0.25*pow((uLstar + uRstar),2) - g*hbar;
            s1s2tilde = fmax(0.0,uLstar*uRstar) - g*hbar;

            // find if sonic problem
            sonic = false;
            if (fabs(s1s2bar) <= criticaltol) sonic = true;
            if (s1s2bar*s1s2tilde <= criticaltol) sonic = true;
            if (s1s2bar*sE1*sE2*E2 <= criticaltol) sonic = true;
            if (fmin(fabs(sE1),fabs(sE2)) < criticaltol) sonic = true;

            // find jump in h, deldelh
            if (sonic)
            {
                deldelh = -delb;
            }
            else
            {
                deldelh = delb*g*hbar/s1s2bar;
            }

            // bounds in case of critical state resonance, or negative states
            if (sE1 < -criticaltol && sE2 > criticaltol)
            {
                deldelh = fmin(deldelh,hstarHLL*(sE2-sE1)/sE2);
                deldelh = fmax(deldelh,hstarHLL*(sE2-sE1)/sE1);
            }
            else if (sE1 >= criticaltol)
            {
                deldelh = fmin(deldelh,hstarHLL*(sE2-sE1)/sE1);
                deldelh = fmax(deldelh, -hL);
            }
            else if (sE2 <= -criticaltol)
            {
                deldelh = fmin(deldelh,hR);
                deldelh = fmax(deldelh, hstarHLL*(sE2-sE1)/sE2);
            }

            // find jump in phi, deldelphi
            if (sonic)
            {
                deldelphi = -g*hbar*delb;
            }
            else
            {
                deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar;
            }

            // bounds in case of critical state resonance, or negative states
            deldelphi = fmin(deldelphi, g*fmax(-hLstar*delb,-hRstar*delb));
            deldelphi = fmax(deldelphi, g*fmin(-hLstar*delb,-hRstar*delb));

            // --- determine fwaves -----------------------------------------------------------

            // first decomposition
            delhdecomp = delh - deldelh;
            delalpha1 = (sE2*delhdecomp - delhu)/(sE2-sE1) - alpha1;
            alpha1 = alpha1 + delalpha1;
            delalpha2 = (delhu - sE1*delhdecomp)/(sE2-sE1) - alpha2;
            alpha2 = alpha2 + delalpha2;

            // second decomposition
            delphidecomp = delphi - deldelphi;
            beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1);
            beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1);

            if ((pow(delalpha2,2) + pow(delalpha1,2)) < pow(convergencetol,2)) break;

            if (sE2 > 0.0 && sE1 < 0.0)
            {
                hLstar = hL + alpha1;
                hRstar = hR - alpha2;
                hustar = huL + beta1;
            }
            else if (sE1 >= 0.0)
            {
                hLstar = hL;
                hustar = huL;
                hRstar = hR - alpha1 - alpha2;
            }
            else if (sE2 <= 0.0)
            {
                hLstar = hR;
                hustar = huR;
                hRstar = hL + alpha1 + alpha2;
            }

            if (hLstar > drytol)
            {
                uLstar = hustar/hLstar;
            }
            else
            {
                hLstar = fmax(hLstar,0.0);
                uLstar = 0.0;
            }

            if (hRstar > drytol)
            {
                uRstar = hustar/hRstar;
            }
            else
            {
                hRstar = fmax(hRstar,0.0);
                uRstar = 0.0;
            }
        }
    }

    delhdecomp = delh - deldelh;
    delphidecomp = delphi - deldelphi;

    // first decomposition
    alpha1 = (sE2*delhdecomp - delhu)/(sE2-sE1);
    alpha2 = (delhu - sE1*delhdecomp)/(sE2-sE1);

    // second decomposition
    beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1);
    beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1);

    //  1st nonlinear wave
    fw[0][0] = alpha1*sE1;
    fw[1][0] = beta1*sE1;
    fw[2][0] = fw[0][0]*vL;

    //  2nd nonlinear wave
    fw[0][2] = alpha2*sE2;
    fw[1][2] = beta2*sE2;
    fw[2][2] = fw[0][2]*vR;

    //  advection of transverse wave
    fw[0][1] = 0.0;
    fw[1][1] = 0.0;
    fw[2][1] = hR*uR*vR - hL*uL*vL - fw[2][1] - fw[2][2];

    //  speeds
    sw[0] = sE1;
    sw[1] = 0.5*(sE1+sE2);
    sw[2] = sE2;
}// === end function Riemann_ssqfwave ========================================================


// === Begin fuction Riemann_fwaves========================================================
// @description: - Solves swe give single left and right states
//               - Solution has two waves
//               - flux - source is decomposed
void riemann_fwaves(int meqn, int mwaves, double hL, double hR, double huL, double huR, 
                    double hvL, double hvR, double bL, double bR, double uL, double uR, 
                    double vL, double vR, double phiL, double phiR, double s1, double s2, double drytol, double g, double sw[], double fw[][])
{
    // Local variables
    double delh, delhu, delphi, delb, delhdecomp, delphidecomp, deldelh, deldelphi;
    double beta1, beta2;

    // determine del vectors
    delh = hR - hL;
    delhu = huR - huL;
    delphi = phiR - phiL;
    delb = bR - bL;

    deldelphi = -g*0.5*(hR + hL)*delb;
    delphidecomp = delphi - deldelphi;

    // flux decomposition
    beta1 = (s2*delhu - delphidecomp)/(s2 - s1);
    beta2 = (delphidecomp - s1*delhu)/(s2 - s1);

    sw[0] = s1;
    sw[1] = 0.5*(s1 + s2);
    sw[2] = s2;

    // 1st nonlinear wave
    fw[0][0] = beta1;
    fw[1][0] = beta1*s1;
    fw[2][0] = beta*vL;

    // 2nd nonlinear wave
    fw[0][2] = beta2;
    fw[1][2] = beta2*s2;
    fw[2][2] = beta2*vR;

    // advsction of transverse wave
    fw[0][1] = 0.0;
    fw[1][1] = 0.0;
    fw[2][1] = hR*uR*vR - hL*uL*vL -fw[2][0] - fw[2][2];
}
// === end function Riemann_fwaves ========================================================


// === Begin fuction Riemann type===========================================================
// @description: Determines the Riemann structure (wave-type in each family)
void riemanntype(double hL, double hR, double uL, double uR, double hm, 
                 double s1m, double s2m, bool rare1, bool rare2, int maxiter,
                 double drytol, double g)
{
    // Local variables
    double u1m, u2m, delu;
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
            for (iter=0; iter <= maxiter; ++iter)
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
            for(iter=0; iter <= maxiter; ++iter)
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
// === end function Riemanntype ==================================================