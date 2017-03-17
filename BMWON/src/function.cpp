#include "initialisationEtDichotomie.h"
#include "function.h"
#include "deriveesEtEquations.h"
#define IterMAX 1000


void CoeffCheb(vector<double> const& fff, vector<double> & coef_fff, chebvar XXX,int OmpOnOff)
{
    double fac = 2.0/XXX.n;

    //On actualise les coefficients chebyshev
    #pragma omp parallel for if(OmpOnOff)
    for(int jjj=0;jjj<XXX.n;jjj++)
    {
        double sum = 0.0;
        for(int kkk=0;kkk<XXX.n;kkk++) sum += fff[kkk]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
        if(jjj!=0) coef_fff[jjj] = fac*sum;
        else coef_fff[jjj] = sum/XXX.n;
    }
}


void CoeffCheb2(vector < vector<double> > const& fff, vector < vector<double> > & coef_fff, int iii, chebvar XXX)
{
    double fac = 2.0/XXX.n;

    //On actualise les coefficients chebyshev

    for(int jjj=0;jjj<XXX.n;jjj++)
    {
        double sum = 0.0;
        for(int kkk=0;kkk<XXX.n;kkk++)
        {
            if(kkk<XXX.n/2) sum += fff[kkk][iii]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
            else sum += fff[XXX.n-kkk-1][iii]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
        }
        if(jjj!=0) coef_fff[iii][jjj] = fac*sum;
        else coef_fff[iii][jjj] = sum/XXX.n;
    }
}

void multiCoeffCheb2_2(vector < vector<double> > const& fff1, vector < vector<double> > & coef_fff1, vector < vector<double> > const& fff2, vector < vector<double> > & coef_fff2, int iii, chebvar XXX)
{
    double fac = 2.0/XXX.n;

    //On actualise les coefficients chebyshev

    for(int jjj=0;jjj<XXX.n;jjj++)
    {
        double sum_1 = 0.0, sum_2 = 0.0;
        for(int kkk=0;kkk<XXX.n;kkk++)
        {
            if(kkk<XXX.n/2)
            {
                sum_1 += fff1[kkk][iii]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
                sum_2 += fff2[kkk][iii]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
            }
            else
            {
                sum_1 += fff1[XXX.n-kkk-1][iii]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
                sum_2 += fff2[XXX.n-kkk-1][iii]*cos(M_PI*jjj*(kkk+0.5)/XXX.n);
            }
        }
        if(jjj!=0)
        {
            coef_fff1[iii][jjj] = fac*sum_1;
            coef_fff2[iii][jjj] = fac*sum_2;
        }
        else
        {
            coef_fff1[iii][jjj] = sum_1/XXX.n;
            coef_fff2[iii][jjj] = sum_2/XXX.n;
        }
    }
}


void CoeffChebgamma(GammaAndDerivatives & gamma,int OmpOnOff)
{
    //On actualise les coefficients chebyshev
    #pragma omp parallel for if(OmpOnOff)
    for(int iii = 0 ; iii <= nrho; iii++) CoeffCheb2(gamma.YYYa_p_rho, gamma.chebcoef_YYYa_rho_p,iii, gamma.Pcheb),CoeffCheb2(gamma.YYYb_p_rho, gamma.chebcoef_YYYb_rho_p,iii, gamma.Pcheb);
}

double chebev(vector<double> const& coef_fff,chebvar const& XXX, double x)
//Chebyshev evaluation: All arguments are input. c[0..m-1] is an array of Chebyshev coefficients, the first m elements of c output from chebft (which must have been called with the same a and b). The Chebyshev polynomia "m-1 k=0 ckTk(y) − c0/2 is evaluated at a point y = [x − (b + a)/2]/[(b− a)/2], and the result is returned as the function value.
{
////    double sum = 0.0;
////    for(int iii = 0; iii < XXX.nev; iii++) sum += coef_fff[iii]*cos(iii*acos((2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)));
//////    sum -= 0.5*coef_fff[0];
////    return sum;

    double dd=0.0,ddd=0.0,sv,y,y2;
    int j;
    y2=2.0*(y=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    for (j=XXX.nev-1;j>=1;j--)// Clenshaw’s recurrence.
    {
        sv=dd;
        dd=y2*dd-ddd+coef_fff[j];
        ddd=sv;
    }
//    return y*dd-ddd+0.5*coef_fff[0];// Last step is different.
    return y*dd-ddd+coef_fff[0];// Last step is different.
}


void multichebev_4(vector<double> const& coef_fff1,vector<double> const& coef_fff2,vector<double> const& coef_fff3,vector<double> const& coef_fff4, double & yy1, double & yy2, double & yy3, double & yy4,chebvar const& XXX, double x)
{//Evaluation in x of 4 functions approximated by chebyshev polynom
    double dd_1=0.0,ddd_1=0.0,sv_1,y_1,y2_1;
    double dd_2=0.0,ddd_2=0.0,sv_2,y_2,y2_2;
    double dd_3=0.0,ddd_3=0.0,sv_3,y_3,y2_3;
    double dd_4=0.0,ddd_4=0.0,sv_4,y_4,y2_4;
    int j;
    y2_1=2.0*(y_1=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    y2_2=2.0*(y_2=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    y2_3=2.0*(y_3=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    y2_4=2.0*(y_4=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    for (j=XXX.nev-1;j>=1;j--)// Clenshaw’s recurrence.
    {
        sv_1=dd_1;
        dd_1=y2_1*dd_1-ddd_1+coef_fff1[j];
        ddd_1=sv_1;

        sv_2=dd_2;
        dd_2=y2_2*dd_2-ddd_2+coef_fff2[j];
        ddd_2=sv_2;

        sv_3=dd_3;
        dd_3=y2_3*dd_3-ddd_3+coef_fff3[j];
        ddd_3=sv_3;

        sv_4=dd_4;
        dd_4=y2_4*dd_4-ddd_4+coef_fff4[j];
        ddd_4=sv_4;
    }
//    return y*dd-ddd+0.5*coef_fff[0];// Last step is different.
    yy1=y_1*dd_1-ddd_1+coef_fff1[0];// Last step is different.
    yy2=y_2*dd_2-ddd_2+coef_fff2[0];// Last step is different.
    yy3=y_3*dd_3-ddd_3+coef_fff3[0];// Last step is different.
    yy4=y_4*dd_4-ddd_4+coef_fff4[0];// Last step is different.
}

void multichebev_2(vector<double> const& coef_fff1,vector<double> const& coef_fff2, double & yy1, double & yy2,chebvar const& XXX, double x)
{//Evaluation in x of 2 functions approximated by chebyshev polynom
    double dd_1=0.0,ddd_1=0.0,sv_1,y_1,y2_1;
    double dd_2=0.0,ddd_2=0.0,sv_2,y_2,y2_2;
    int j;
    y2_1=2.0*(y_1=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    y2_2=2.0*(y_2=(2.0*x-XXX.xmin-XXX.xmax)/(XXX.xmax-XXX.xmin)); //Change of variable.
    for (j=XXX.nev-1;j>=1;j--)// Clenshaw’s recurrence.
    {
        sv_1=dd_1;
        dd_1=y2_1*dd_1-ddd_1+coef_fff1[j];
        ddd_1=sv_1;

        sv_2=dd_2;
        dd_2=y2_2*dd_2-ddd_2+coef_fff2[j];
        ddd_2=sv_2;
    }
//    return y*dd-ddd+0.5*coef_fff[0];// Last step is different.
    yy1=y_1*dd_1-ddd_1+coef_fff1[0];// Last step is different.
    yy2=y_2*dd_2-ddd_2+coef_fff2[0];// Last step is different.
}

void Derivation(GammaAndDerivatives& gamma, int iii)
{
        multiCoeffCheb2_2(gamma.YYYa_p_rho, gamma.chebcoef_YYYa_rho_p,gamma.YYYb_p_rho, gamma.chebcoef_YYYb_rho_p,iii, gamma.Pcheb);
//        CoeffCheb2(gamma.YYYa_p_rho, gamma.chebcoef_YYYa_rho_p,iii, gamma.Pcheb);
//        CoeffCheb2(gamma.YYYb_p_rho, gamma.chebcoef_YYYb_rho_p,iii, gamma.Pcheb);
        mulitchder_2(gamma.Pcheb,gamma.chebcoef_YYYa_rho_p[iii],gamma.chebcoef_YYYap1_rho_p[iii],gamma.chebcoef_YYYb_rho_p[iii],gamma.chebcoef_YYYbp1_rho_p[iii]);
//        chder(gamma.Pcheb,gamma.chebcoef_YYYa_rho_p[iii],gamma.chebcoef_YYYap1_rho_p[iii]);
//        chder(gamma.Pcheb,gamma.chebcoef_YYYb_rho_p[iii],gamma.chebcoef_YYYbp1_rho_p[iii]);
        gamma.WWW1_rho[iii] = derivePrm(gamma.WWW_rho,iii,nrho);
        gamma.WWW2_rho[iii] = deriveScd(gamma.WWW_rho,iii,nrho);

        for(int ipp = 0; ipp < nppp; ipp++)
        {
            gamma.JJJ3LL_p_rho[ipp][iii] = 0.0,gamma.JJJ3TT_p_rho[ipp][iii] = 0.0,gamma.JJJ3LT_p_rho[ipp][iii] = 0.0,gamma.JJJ3TL_p_rho[ipp][iii] = 0.0;
            gamma.YYYa1_p_rho[ipp][iii] = derivePrm(gamma.YYYa_p_rho[ipp],iii,nrho);
            gamma.YYYa2_p_rho[ipp][iii] = deriveScd(gamma.YYYa_p_rho[ipp],iii,nrho);
            gamma.YYYb1_p_rho[ipp][iii] = derivePrm(gamma.YYYb_p_rho[ipp],iii,nrho);
            gamma.YYYb2_p_rho[ipp][iii] = deriveScd(gamma.YYYb_p_rho[ipp],iii,nrho);
        }
        multiCoeffCheb2_2(gamma.YYYa1_p_rho, gamma.chebcoef_YYYa1_rho_p,gamma.YYYb1_p_rho, gamma.chebcoef_YYYb1_rho_p,iii, gamma.Pcheb);
//        CoeffCheb2(gamma.YYYa1_p_rho, gamma.chebcoef_YYYa1_rho_p,iii, gamma.Pcheb);
//        CoeffCheb2(gamma.YYYb1_p_rho, gamma.chebcoef_YYYb1_rho_p,iii, gamma.Pcheb);

}//fin fonction

void Integrals(GaussLegendre gausslegvar, GammaAndDerivatives & gamma, int iii)
{
    //Calcul des intégrales
    double bpaq1 = 0, bmaq1 = qInt, bpaq2 = 0.5*qInt, bmaq2 = 0.5*qInt;
    double q1, q2, p, qsqr, facReg, propgl, propgt, auxl, auxt, yya, yyb, yyb1, yya1;


    double mk2 = gamma.WWW_rho[iii] + 2*iii*gamma.WWW1_rho[iii], www = gamma.WWW_rho[iii], www1 = gamma.WWW1_rho[iii];
    double int2ll = 0.0, int2tt = 0.0, int3ll = 0.0, int3tt = 0.0, int3lt = 0.0, int3tl = 0.0, inta = 0.0, int1p = 0.0;
    double yamax, ybmax,ya_p,yb_p;
    multichebev_2(gamma.chebcoef_YYYa_rho_p[iii],gamma.chebcoef_YYYb_rho_p[iii],yamax,ybmax,gamma.Pcheb,p_max);

    for(int jjj=1;jjj<=gausslegvar.ngl;jjj++)
    {
        q2 = (bpaq2+bmaq2*gausslegvar.x[jjj]);
        for(int kkk=1;kkk<=gausslegvar.ngl;kkk++)
        {
            q1 = (bpaq1+bmaq1*gausslegvar.x[kkk]);
            qsqr = pow(q1,2)+pow(q2,2);

            facReg = gausslegvar.w[jjj]*gausslegvar.w[kkk]*pow(q2,gausslegvar.d-2)*sss(qsqr,gausslegvar.alpha,gamma.eta);

            multichebev_4(gamma.chebcoef_YYYa_rho_p[iii],gamma.chebcoef_YYYb_rho_p[iii],gamma.chebcoef_YYYa1_rho_p[iii],gamma.chebcoef_YYYb1_rho_p[iii],yya,yyb,yya1,yyb1,gamma.Pcheb,sqrt(qsqr));

//            yyb = chebev(gamma.chebcoef_YYYb_rho_p[iii],gamma.Pcheb,sqrt(qsqr));
//            yya1 = chebev(gamma.chebcoef_YYYa1_rho_p[iii],gamma.Pcheb,sqrt(qsqr));
//
//            propgl = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*( 1 + chebev(gamma.chebcoef_YYYa_rho_p[iii],gamma.Pcheb,sqrt(qsqr)) + 2*iii*drho*yyb ) + mk2,-1 ) ;
//            propgt = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*( 1 + chebev(gamma.chebcoef_YYYa_rho_p[iii],gamma.Pcheb,sqrt(qsqr)) ) + www,-1 ) ;

            propgl = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*( 1 + yya + 2*iii*drho*yyb ) + mk2,-1 ) ;
            propgt = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*( 1 + yya ) + www,-1 ) ;

            auxl = facReg*pow(propgl,2);
            auxt = facReg*pow(propgt,2);

//            int1p += -(gamma.n-1)*auxt*(qsqr*yya1 + www1)/drho - auxl*(qsqr*( yya1/drho + 2*yyb + 2*iii*chebev(gamma.chebcoef_YYYb1_rho_p[iii],gamma.Pcheb,sqrt(qsqr)) ) + (3*www1 + 2*iii*gamma.WWW2_rho[iii])/drho);
            int1p += -(gamma.n-1)*auxt*(qsqr*yya1 + www1)/drho - auxl*(qsqr*( yya1/drho + 2*yyb + 2*iii*yyb1 ) + (3*www1 + 2*iii*gamma.WWW2_rho[iii])/drho);
            inta += facReg*propgl*propgt*(propgl+propgt)*(qsqr*yyb + www1/drho);
            int2ll += auxl;
            int2tt += auxt;
            int3ll += auxl * propgl;
            int3tt += auxt * propgt;
            int3lt += auxl * propgt;
            int3tl += auxt * propgl;

            for(int ipp=0;ipp<nppp;ipp++)
            {
                p = gamma.Pcheb.xcheb[ipp];
                qsqr = pow(q1,2) + pow(q2,2) + pow(p,2) + 2*p*q1;
                p = sqrt(qsqr);

                if(p<p_max)
                {
                    multichebev_2(gamma.chebcoef_YYYa_rho_p[iii],gamma.chebcoef_YYYb_rho_p[iii],ya_p,yb_p,gamma.Pcheb,p);
                    propgl = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*(1 + ya_p + 2*iii*drho*yb_p) + mk2,-1);
                    propgt = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*(1 + ya_p) + www,-1);
                }
                else propgl = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*( 1 + yamax + 2*iii*drho*ybmax ) + mk2,-1), propgt = pow( regulator(qsqr,gausslegvar.alpha) + qsqr*(1 + yamax) + www,-1);
                gamma.JJJ3LL_p_rho[ipp][iii] += gausslegvar.norme * auxl * propgl;
                gamma.JJJ3TT_p_rho[ipp][iii] += gausslegvar.norme * auxt * propgt;
                gamma.JJJ3LT_p_rho[ipp][iii] += gausslegvar.norme * auxl * propgt;
                gamma.JJJ3TL_p_rho[ipp][iii] += gausslegvar.norme * auxt * propgl;
            }
        }
    }
    gamma.WWWI1_rho[iii] = int1p * gausslegvar.norme;
    gamma.IIIA_rho[iii] = inta * gausslegvar.norme;
    gamma.III2LL_rho[iii] = int2ll * gausslegvar.norme;
    gamma.III2TT_rho[iii] = int2tt * gausslegvar.norme;
    gamma.III3LL_rho[iii] = int3ll * gausslegvar.norme;
    gamma.III3TT_rho[iii] = int3tt * gausslegvar.norme;
    gamma.III3LT_rho[iii] = int3lt * gausslegvar.norme;
    gamma.III3TL_rho[iii] = int3tl * gausslegvar.norme;
}

void DeriveesEtIntegrales(GammaAndDerivatives & gamma, GaussLegendre gausslegvar)
{
    #pragma omp parallel for if(gamma.ompOn)
    for(int iii = 0; iii <= nrho ; iii++) Derivation(gamma,iii), Integrals(gausslegvar,gamma,iii);
}


void flowstepfoward(GammaAndDerivatives & gamma, ofstream & log)
{
    gamma.eta = etastepforward(gamma);
    
    #pragma omp parallel for if(gamma.ompOn)
    for(int iii = 0; iii <= nrho ; iii++)
    {
        gamma.WWW_rho[iii] += dtt*flowOfW(gamma,iii);
        for(int ipp = 0; ipp < nppp; ipp++)
        {
            gamma.YYYa_p_rho[ipp][iii] += dtt*flowOfYa(gamma,ipp,iii);
            gamma.YYYb_p_rho[ipp][iii] += dtt*flowOfYb(gamma,ipp,iii);
        }
    }
    
    gamma.Zk -= gamma.eta*gamma.Zk*dtt;
    gamma.YYYa_p_rho[gamma.pPrescription][gamma.rhoPrescription]=0.0;
}


double regulator(double y, double alpha)
{
    return alpha*y/(exp(y)-1);
}


double sss(double y, double alpha, double eta)
{
    return alpha*y*(-eta+2*y*exp(y)/(exp(y)-1))/(exp(y)-1);
}

