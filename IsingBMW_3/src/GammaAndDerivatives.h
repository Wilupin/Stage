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


/* --------------------------------------------------------
   Description de la classe : 
   Cette classe est la classe centrale du code IsingBMW. 
   Elle permet de calculer la decomposition de Chebychev, 
   les derivees des fonctions et les integrales a l'aide 
   de la fonction ChebychevDeriveesEtIntegrales. Ensuite,
   une fois tout ces calculs realise elle permet de calculer
   les equations de flots et donc de remplacer les fonctions 
   par leur nouvel itere temporel en appelant FlowStepForward.
   -------------------------------------------------------- */


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
  double const& get_tt(){return tt;};
  int const& get_ompOn(){return ompOn;};
  int const& get_itMax() {return itMax;}; 
  int const& get_irho0(){return irho0;};
  int const& get_nrho() {return nrho;};
  double const& get_drho() {return drho;};
  double const& get_dt() {return dtt;}; 
  
  
  // Setteur
  void set_kktt(int const& n_t);

  void set_V_rho(double const& val, int ir) {V_rho[ir] = val;};
  void set_W_rho(double const& val, int ir) {W_rho[ir] = val;};
  void set_Delta_px_py_rho(double const& val, int ipx, int ipy, int ir)
  {Delta_px_py_rho[ipx][ipy][ir] = val;};
  
  
  // Fonctions pour les sorties;
  void set_irho0();
  void SortieDiversesGnuPlot(int nrun);
  void SortieDataRho0();
  void SortieLogRho0();

  // Equations de flot
  double FlowOfDeltaDim(int ipx, int ipy, int ir);
  double FlowOfWDim(int ir);
  double FlowOfVDim(int ir);

  // Pas de temps euler, avancee du flot
  void FlowStepForward2D();

  // Integrale et derivation
  void ChebychevDeriveesEtIntegrales();
  void DerivationDim2D(int ir);
  void IntegralsDim2D(int ir);
  void UpdatePropagatorQ(double qx, double qy, int ir);
  double PropagatorPQ(double pxqx, double pyqy, int ir);
  double fToIntI(double qx, double qy, int ir, int i);
  double fToIntI11(double qx, double qy, int ir);
  double fToIntJ(double qx, double qy, double pxpqx, double pypqy,
		 double pxmqx, double pymqy, int ir, int i);
  double RegulatorDim(double qx, double qy);
  double DerRegulatorDim(double qx, double qy);
  double e0(double qx, double qy);
  void Recenter(double px, double py, double qx, double qy,
		double & pxpqx, double & pypqy, double & pxmqx, double & pymqy);
  void DerivationI11Dim2D(int ir);

  
  
 private:
  
  Chebychev Cheb;
  Derivateur Der;
  GaussLegendre GL;

  // Pointeur vers le gestionnaire des entrees et sorties
  InOut *ptrInOut = InOut::getInstance();
  
  // Variables for the resolution
  double kk, tt, dtt;
  int itt, itMax, tMax, tMin;
  
  double d, drho, qMax, alpha, pMax, mu;
  int ompOn, irho0, nrho;

  // Variables de stockage
  std::vector<double> propagatorQ;
  
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
