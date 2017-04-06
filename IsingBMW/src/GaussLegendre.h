#ifndef DEF_GAUSSLEGENDRE
#define DEF_GAUSSLEGENDRE

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>



#include "Parametres.h"

class GaussLegendre
{
  
 public:
  
  // Constructeur et destructeur
  GaussLegendre(int nd);
  ~GaussLegendre(){};

  // Accesseur
  int const& get_ngl(){return ngl;};
  double const& get_norme(){return norme;};
  double const& get_w(int iii){return w[iii];};
  double const& get_x(int iii){return x[iii];};
  
  // Fonction de classe
  double normalisationIntegrale(double const& d);

 private:
  int ngl;
  double d, norme, qMax;
  double EPS = 3.0e-14;
  std::vector<double> w,x;

};


#endif
