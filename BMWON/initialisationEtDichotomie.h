#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>

#define p_max 6
#define p_min -6
#define nppp 20
#define npcheb 2*nppp
#define npchebev int(3*npcheb/4)

#define rho_min 0.0
#define rho_max 10.0
#define drho 0.1
#define nrho int((rho_max-rho_min)/drho)


#define dtt -1.0e-4
#define pprecision 1e-15
#define qInt 4
#define ngausslegendre 40
#define EPS 3.0e-14
#define trigavc 1e-10


using namespace std;

typedef struct DichoParameter
{
    double r0, u0, eta_avc, w0_avc, Zk_avc;
    double rmin, rmax, rht, rlt, seuil;
    int OnOff, tht, tlt, nrun, avc, tttavc;
    vector < vector<double> > YYYa_p_rho_avc,YYYb_p_rho_avc;
    vector<double> WWW_rho_avc;
}DichoParameter;

typedef struct chebvar
{
    int n, nev;
    double xmin, xmax;
    vector<double> xcheb;
}chebvar;


typedef struct GammaAndDerivatives
{
    vector < vector<double> > YYYa_p_rho, YYYa1_p_rho,YYYa2_p_rho,chebcoef_YYYa1_rho_p,chebcoef_YYYa_rho_p, chebcoef_YYYap1_rho_p;
    vector < vector<double> > YYYb_p_rho, YYYb1_p_rho,YYYb2_p_rho,chebcoef_YYYb1_rho_p,chebcoef_YYYb_rho_p, chebcoef_YYYbp1_rho_p;
    vector < vector<double> > JJJ3LL_p_rho, JJJ3TT_p_rho, JJJ3LT_p_rho, JJJ3TL_p_rho;
    vector<double> III2LL_rho, III2TT_rho, III3LL_rho, III3TT_rho, III3LT_rho, III3TL_rho, IIIA_rho;
    vector<double> WWW_rho, WWW1_rho, WWW2_rho, WWWI1_rho;
    chebvar Pcheb;
    double d, eta, Zk, n, kk;
    int pPrescription,rhoPrescription, tmax,ompOn, irho0;
}GammaAndDerivatives;

typedef struct GaussLegendre
{
    int ngl;
    double d, alpha, norme;
    vector<double> w, x;
}GaussLegendre;


void ChebRoot(chebvar& XXX, char OutFile[], int ncheb, int nevcheb, double x_min, double x_max);

void gaulegWeightAbscissas(GaussLegendre& gausslegvar);
double normalisationIntegrale(double d);
void initialiseGamma(GammaAndDerivatives& gamma);

void ReadAndInitialise(GammaAndDerivatives & gamma, GaussLegendre & gausslegvar, DichoParameter & dichoprm);

int initialising(GammaAndDerivatives& gamma, DichoParameter & dichoprm);
int dichotomie(double w0, DichoParameter & dichoprm, int ttt, ofstream & log);
int activationDichoSup(GammaAndDerivatives const& gamma,DichoParameter & dichoprm, ofstream & log);

void sortieDiversesMathematica(GammaAndDerivatives const& gamma);
void sortieDiversesGnuPlot(GammaAndDerivatives const& gamma, int ttt, int nrun);

inline int set_irho_0(vector<double> www)
{
    double wi=www[0],wip1=www[1];
    int iii = 0;
    while(!(wi<=0&&wip1>0)&&iii<nrho) iii++, wi=www[iii], wip1=www[iii+1];
    return iii;
}

inline void OutEtaW0(int ttt, double eta, double w0, int irho0)
{
    ofstream w0etaout("data/w0eta",ios::app);
    w0etaout.precision(16);
    w0etaout<<ttt*dtt<<"\t"<<eta<<"\t"<<w0<<"\t"<<irho0<<endl;
    w0etaout.close();
}

inline double SetKd(double d)
{
    return pow(2,1-d)*pow(M_PI,-d/2)*pow(d*tgamma(d/2),-1);
}

