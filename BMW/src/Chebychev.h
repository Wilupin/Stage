#ifndef DEF_CHEBYCHEV
#define DEF_CHEBYCHEV

/* ------------------------------
** 
** Gaetan Facchinetti
** Mail : gaetanfacc@gmail.com
**
** Created : 21/03/2017
**
** -----------------------------*/


/* Chebychev : Une classe SINGLETON qui defini les proprietes
** de la decomposition sur les polynomes de Chebychev et 
** l'effectue a l'aide de differentes fonctions
*/


#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>


class Chebychev
{

 public:

  // Fonctions qui en font un singleton
  static Chebychev *getInstance()
  {
    if(NULL == _singletonChebychev)
      {
	_singletonChebychev = new Chebychev;
      }
    return _singletonChebychev;
  }

  static void kill()
  {
    if(NULL != _singletonChebychev)
      {
	delete _singletonChebychev;
	_singletonChebychev = NULL;
      }
  }

  
  // Fonction de classe
  void   ChebRoot(std::string OutFile, int ncheb, int nevcheb, double x_min, double x_max);
  void   CoeffCheb(std::vector<double> const& fff,
		   std::vector<double> & coef_fff, int OmpOnOff);
  void   CoeffCheb2(std::vector < std::vector<double> > const& fff,
		    std::vector < std::vector<double> > & coef_fff, int iii);
  void   multiCoeffCheb2_2(std::vector < std::vector<double> > const& fff1,
			   std::vector < std::vector<double> > & coef_fff1,
			   std::vector < std::vector<double> > const& fff2,
			   std::vector < std::vector<double> > & coef_fff2, int iii);
  double chebev(std::vector<double> const& coef_fff, double const& x) const;
  void   multichebev_4(std::vector<double> const& coef_fff1,
		       std::vector<double> const& coef_fff2,
		       std::vector<double> const& coef_fff3,
		       std::vector<double> const& coef_fff4,
		       double & yy1, double & yy2, double & yy3,
		       double & yy4, double x);
  void multichebev_2(std::vector<double> const& coef_fff1,
		     std::vector<double> const& coef_fff2,
		     double & yy1, double & yy2, double x);

  
  // Accesseur
  int const& get_n() const {return n;}
  double const& get_xmin() const {return xmin;}
  double const& get_xmax() const {return xmax;}
  double const& get_xcheb(int i) const {return xcheb[i];}

  // Setteur
  void set_n(int const& n_n){n = n_n;}
  void set_xmin(double const& n_xmin){xmin = n_xmin;}
  void set_xmax(double const& n_xmax){xmax = n_xmax;}
  void set_xcheb(std::vector<double> const& n_xcheb){xcheb = n_xcheb;} 
  
  

 private:

  // Constructeur et destructeur sont privees (c'est un singleton)
  Chebychev(){};
  ~Chebychev(){};
  
  static Chebychev *_singletonChebychev;
  int n, nev;
  double xmin, xmax;
  std::vector<double> xcheb;
  
};


#endif
