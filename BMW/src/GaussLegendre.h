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



class GaussLegendre
{
  
 public:
  
  // Constructeur et destructeur
  GaussLegendre(){ngl = ngausslegendre;};
  ~GaussLegendre(){};

  // Accesseur
  int const& get_ngl(){return ngl;};
  double const& get_d(){return d;};
  double const& get_alpha(){return alpha;};
  double const& get_norme(){return norme;};
  double const& get_w(int iii){return w[iii];};
  double const& get_x(int iii){return x[iii];};

  // Setteur
  void set_alpha(double const& n_alpha){alpha = n_alpha;};
  void set_d(double const& n_d)        {d = n_d;};
  
  // Fonction de classe
  static void gaulegWeightAbscissas(GaussLegendre & gausslegvar);
  static double normalisationIntegrale(double const& d);
  static double SetKd(double const& d);

 private:
  int ngl;
  double d, alpha, norme;
  std::vector<double> w,x;

};


#endif
