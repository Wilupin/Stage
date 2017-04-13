#ifndef DEF_DERIVATEUR
#define DEF_DERIVATEUR


#include "Chebychev.h"


class Derivateur
{

 public:
  
  // Constructeur et destructeurs
  Derivateur(){};
  Derivateur(Chebychev* nCheb){Cheb = nCheb;};
  ~Derivateur(){};
  
  
  // Fonction de classe
  double der1bulk(double, double, double,double, double);
  double der1FstPt(double, double, double,double, double);
  double der1ScdPt(double, double, double,double, double);
  double der1LastPtm1(double, double, double,double, double);
  double der1LastPt(double, double, double,double, double);

  double der2bulk(double, double, double,double, double);
  double der2FstPt(double, double, double,double, double);
  double der2ScdPt(double, double, double,double, double);
  double der2LastPtm1(double, double, double,double, double);
  double der2LastPt(double, double, double,double, double);


  void chder(std::vector<double> const& c, std::vector<double> & cder);
  void multichder_2(std::vector<double> const& c1, std::vector<double> & cder1,
		     std::vector<double> const& c2, std::vector<double> & cder2);
  void chder2(std::vector < std::vector<double> > const& c,
	      std::vector < std::vector<double> > & cder, int iii);
  double derivePrm(std::vector<double> const& fff, int iii, int nbPtRho);
  double deriveScd(std::vector<double> const& fff, int iii, int nbPtRho);

  
 private:
  Chebychev *Cheb; 

  
};

#endif
