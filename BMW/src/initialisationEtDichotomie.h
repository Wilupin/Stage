#include "GammaAndDerivatives.h"

using namespace std;


// Structure :: data for dichotomie procedure
typedef struct DichoParameter
{
  double r0, u0, eta_avc, w0_avc, Zk_avc;
  double rmin, rmax, rht, rlt, seuil;
  int OnOff, tht, tlt, nrun, avc, tttavc;
  vector < vector<double> > YYYa_p_rho_avc,YYYb_p_rho_avc;
  vector<double> WWW_rho_avc;
}DichoParameter;



void ReadAndInitialise(GammaAndDerivatives & gamma,
		       GaussLegendre & gausslegvar,
		       DichoParameter & dichoprm);

int initialising(GammaAndDerivatives& gamma, DichoParameter & dichoprm);
int dichotomie(double w0, DichoParameter & dichoprm, int ttt, ofstream & log);


inline int set_irho_0(vector<double> www)
{
    double wi=www[0],wip1=www[1];
    int iii = 0;
    while(!(wi<=0&&wip1>0)&&iii<nrho) iii++, wi=www[iii], wip1=www[iii+1];
    return iii;
}

inline void OutEtaW0(int ttt, double eta, double w0, int irho0)
{
    ofstream w0etaout("../data/w0eta",ios::app);
    w0etaout.precision(16);
    w0etaout<<ttt*dtt<<"\t"<<eta<<"\t"<<w0<<"\t"<<irho0<<endl;
    w0etaout.close();
}



