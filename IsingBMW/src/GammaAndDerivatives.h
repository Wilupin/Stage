#ifndef DEF_GAMMAANDDERIVATIVES
#define DEF_GAMMAANDDERIVATIVES

/* ------------------------------
** 
** Gaetan Facchinetti
** Mail : gaetanfacc@gmail.com
**
** Created : 21/03/2017
**
** -----------------------------*/


#include "Chebychev.h"
#include "GaussLegendre.h"
#include "Derivateur.h"
#include "InOut.h"


class GammaAndDerivatives
{
 public:
  
  // Constructeur et destructeur
  GammaAndDerivatives();
  ~GammaAndDerivatives(){};


  // Accesseur
  double const& get_d(){return d;};
  double const& get_kk(){return kk;};
  int const& get_ompOn(){return ompOn;};
  int const& get_tmax() {return tmax;}; 
  int const& get_irho0(){return irho0;};
  int const& get_nrho() {return nrho;};
  double const& get_drho() {return drho;};
  double const& get_dt() {return dtt;}; 
  
  
  // Setteur
  void set_d(double const& n_d)           {d = n_d;};
  void set_kk(double const& n_k)          {kk = n_k;};
  void set_ompOn(int const& n_o)          {ompOn = n_o;};
  void set_tmax(int const& n_t)           {tmax = n_t;};
  
  // Fonctions autres;
  void set_irho0();
  void sortieDiversesGnuPlot(int ttt, int nrun);

  
  /* -------------------------------------------------
     Fonctions pour le calcul dimensionne Ising d = 2 
     ------------------------------------------------- */

  // Equations de flot
  double flowOfDeltaDim(int ipx, int ipy, int ir);
  double flowOfWDim(int ir);
  double flowOfVDim(int ir);

  // Calcul des coefficients de chebichev
  void CoeffChebGamma2D();

  // Pas de temps euler
  void flowStepForward2D(std::ofstream & log);

  // Integrale et derivation
  void DeriveesEtIntegrales();
  void DerivationDim2D(int ir);
  void IntegralsDim2D(int ir);
  void updatePropagatorQ(double qx, double qy, int ir);
  void updatePropagatorPQ(double pxpqx, double pypqy, int ir);
  double fToIntI(double qx, double qy, int i);
  double fToIntJ(double qx, double qy, int i);
  double RegulatorDim(double qx, double qy);
  double DerRegulatorDim(double qx, double qy);
  double e0(double qx, double qy);
  void Recenter(double px, double py, double qx, double qy,
		double & npxpqx, double & npypqy);
  void DerivationI11Dim2D(int ir);

  
 private:
  
  Chebychev Cheb;
  Derivateur Der;
  GaussLegendre GL;

  // Pointeur vers le gestionnaire des entrees et sorties
  InOut *ptr = InOut::getInstance();
  
  // Variables for the resolution
  double d, kk, drho, qMax, dtt, alpha, pMax, mu;
  int tmax, ompOn, irho0, nrho;

  
  /* -------------------------------------------------
     Variables pour le calcul dimensionne Ising d = 2 
     ------------------------------------------------- */

  // Variables de stockage
  double uk, gammaq, ek, m2k, valDeltaQ, valDeltaPQ;
  double propagatorQ, propagatorPQ;
  
  // Nouvelles fonctions (Variables dimensionnees)
  std::vector <std::vector < std::vector <double> > > Delta_px_py_rho, ChebDelta_px_py_rho;
  std::vector <std::vector < std::vector <double> > > Delta1_px_py_rho, ChebDelta1_px_py_rho;
  std::vector <std::vector < std::vector <double> > > Delta2_px_py_rho;
  std::vector <double> W_rho;            // W_rho = dVk(rho)/drho
  std::vector <double> V_rho;
  std::vector <double> W1_rho, W2_rho;   // Derivees

  // Integrales
  std::vector< std::vector < std::vector<double>>> J3_px_py_rho;
  std::vector<double> I1_rho, I2_rho, I3_rho;
  std::vector<double> I11_rho;
};





#endif
