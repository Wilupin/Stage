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



class GammaAndDerivatives
{
 public:
  
  // Constructeur et destructeur
  GammaAndDerivatives(GaussLegendre* m_gausslegvar);
  ~GammaAndDerivatives(){delete gausslegvar, XXX;};


  // Accesseur
  double const& get_d(){return d;};
  double const& get_n(){return n;};
  double const& get_eta(){return eta;};
  double const& get_Zk(){return Zk;};
  double const& get_kk(){return kk;};
  int const& get_ompOn(){return ompOn;};
  int const& get_tmax() {return tmax;}; 
  int const& get_pPrescription(){return pPrescription;};
  int const& get_rhoPrescription(){return rhoPrescription;};
  int const& get_irho0(){return irho0;};

  double const& get_WWW_rho(int iii){return WWW_rho[iii];};
  std::vector<double> const& get_WWW_rho_vec(){return WWW_rho;};


  
  // Setteur
  void set_d(double const& n_d)           {d = n_d;};
  void set_n(double const& n_n)           {n = n_n;};
  void set_eta(double const& n_e)         {eta = n_e;};
  void set_Zk(double const& n_z)          {Zk = n_z;};
  void set_kk(double const& n_k)          {kk = n_k;};
  void set_ompOn(int const& n_o)          {ompOn = n_o;};
  void set_tmax(int const& n_t)           {tmax = n_t;};
  void set_pPrescription(int const& n_p)  {pPrescription = n_p;};
  void set_rhoPrescription(int const& n_r){rhoPrescription = n_r;};
  void set_irho0(int const& n_r)           {irho0 = n_r;};

  void set_WWW_rho(double const& n_W, int iii) {WWW_rho[iii] = n_W;};
  void set_YYYa_p_rho(double const& n_y, int ipp, int  iii) {YYYa_p_rho[ipp][iii] = n_y;};
  void set_YYYb_p_rho(double const& n_y, int ipp, int  iii) {YYYb_p_rho[ipp][iii] = n_y;};



  
  // Fonctions de classe
  void flowstepfoward(std::ofstream & log);
  double flowOfYa(int ipp, int iii);
  double flowOfYb(int ipp, int iii);
  double flowOfW(int iii);
  double etastepforward();

  void Derivation(int iii);
  void Integrals(int iii);
  void DeriveesEtIntegrales();
  void CoeffChebgamma();

  static double regulator(double y, double alpha);
  static double sss(double y, double alpha, double eta);
    
  void sortieDiversesGnuPlot(int ttt, int nrun);
  
    
  
    
 private:

  Chebychev *XXX             = Chebychev::getInstance();
  Derivateur *Der            = Derivateur::getInstance();
  GaussLegendre *gausslegvar = new GaussLegendre();
  
  // Two variables dependant functions of p and rho
  std::vector < std::vector<double> > YYYa_p_rho, YYYa1_p_rho,YYYa2_p_rho;
  std::vector < std::vector<double> > chebcoef_YYYa1_rho_p,chebcoef_YYYa_rho_p, chebcoef_YYYap1_rho_p;
  std::vector < std::vector<double> > YYYb_p_rho, YYYb1_p_rho,YYYb2_p_rho;
  std::vector < std::vector<double> > chebcoef_YYYb1_rho_p,chebcoef_YYYb_rho_p, chebcoef_YYYbp1_rho_p;
  std::vector < std::vector<double> > JJJ3LL_p_rho, JJJ3TT_p_rho, JJJ3LT_p_rho, JJJ3TL_p_rho;

  // One variable dependant functions of rho
  std::vector<double> III2LL_rho, III2TT_rho, III3LL_rho;
  std::vector<double> III3TT_rho, III3LT_rho, III3TL_rho, IIIA_rho;
  std::vector<double> WWW_rho, WWW1_rho, WWW2_rho, WWWI1_rho;

  // Other variables
  double d, eta, Zk, n, kk;
  int pPrescription,rhoPrescription, tmax, ompOn, irho0;

};





#endif
