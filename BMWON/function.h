
// Fonctions Chebyshev
void CoeffChebgamma(GammaAndDerivatives & gamma,int OmpOnOff);
void CoeffCheb(vector<double> const& fff, vector<double> & coef_fff, chebvar XXX,int OmpOnOff);
void CoeffCheb2(vector < vector<double> > const& fff, vector < vector<double> > & coef_fff, int iii, chebvar XXX);
void multiCoeffCheb2_2(vector < vector<double> > const& fff1, vector < vector<double> > & coef_fff1, vector < vector<double> > const& fff2, vector < vector<double> > & coef_fff2, int iii, chebvar XXX);
double chebev(vector<double> const& coef_fff,chebvar const& XXX, double x);
void multichebev_4(vector<double> const& coef_fff1,vector<double> const& coef_fff2,vector<double> const& coef_fff3,vector<double> const& coef_fff4, double & y1, double & y2, double & y3, double & y4,chebvar const& XXX, double x);
void multichebev_2(vector<double> const& coef_fff1,vector<double> const& coef_fff2, double & yy1, double & yy2,chebvar const& XXX, double x);

//Fonctions pour intégrales gauss-legendre
double regulator(double y, double alpha);
double sss(double y, double alpha, double eta);


void Derivation(GammaAndDerivatives& gamma, int iii);
void Integrals(GaussLegendre gausslegvar, GammaAndDerivatives & gamma, int iii);
void DeriveesEtIntegrales(GammaAndDerivatives & gamma, GaussLegendre gausslegvar);
void flowstepfoward(GammaAndDerivatives & gamma, ofstream & log);



