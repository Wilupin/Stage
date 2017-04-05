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
  double const& get_alpha(){return alpha;};
  double const& get_norme(){return norme;};
  double const& get_w(int iii){return w[iii];};
  double const& get_x(int iii){return x[iii];};

  int const& get_intTriangleY(int iii){return intTriangleY[iii];};
  
  
  // Fonction de classe
  double normalisationIntegrale(double const& d);

 private:
  int ngl;
  double d, alpha, norme, qMax;
  double EPS = 3.0e-14;
  std::vector<double> w,x;

  /* Vecteur qui contient les points d'integration que l'on doit regarder
     lorsque l'on fait une integrale sur un triangle (par symmetrie dans le
     carre) afin de ne pas avoir a compter deux fois les mêmes points) */
  std::vector<int> intTriangleY;

};


#endif
