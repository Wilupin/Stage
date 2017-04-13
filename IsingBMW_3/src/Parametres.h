#ifndef DEF_PARAMETRES
#define DEF_PARAMETRES


#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

class Parametres
{
    
public:

  // Fonctions qui en font un singleton
  static Parametres *getInstance()
  {
    if(NULL == _Parametres)
      {
	_Parametres = new Parametres;
      }
    return _Parametres;
  }

  static void kill()
  {
    if(NULL != _Parametres)
      {
	delete _Parametres;
	_Parametres = NULL; 
      }
  }

  
  // Accesseurs
  double const& get_rhoMin() const {return rhoMin;};
  double const& get_rhoMax() const {return rhoMax;};
  double const& get_drho()   const {return drho;};
  double const& get_qMax()   const {return qMax;};
  double const& get_pMax()   const {return pMax;};
  double const& get_pMin()   const {return pMin;};
  
  int const& get_npGL() const {return npGL;};
  int const& get_nrho() const {return nrho;};

  int const& get_npCheb()   const {return npCheb;};
  int const& get_npChebev() const {return npChebev;};
  
  double const& get_alpha() const {return alpha;};
  double const& get_mu()    const {return mu;};

  double const& get_beta()    const {return beta;};
  
  
  double const& get_tMax()  const {return tMax;};
  double const& get_tMin()  const {return tMin;};
  double const& get_dt()    const {return dt;};

  int const& get_tMaxInt()  const {return tMaxInt;};
    
  int const& get_ompOn()           const {return ompOn;};
  int const& get_numberOfThreads() const {return numberOfThreads;};

  int const& get_itSaveInLog()       const {return itSaveInLog;}
  int const& get_itSaveInData()      const {return itSaveInData;}
  int const& get_itSaveEvthgInData() const {return itSaveEvthgInData;}
  
  
  // Fonctions de classe 
  std::string trim(std::string const& str);


    
 private:

   // Constructeur et destructeur
  Parametres();
  ~Parametres(){};

  static Parametres *_Parametres;

  
  // Parametres de calcul 
  double rhoMin, rhoMax, drho, qMax;
  int npGL, nrho;
  
  // Parametres du probleme
  double alpha, mu;

  // Parametre temperature
  double beta;
  
  // Parametres du temps
  double tMax, tMin, dt;
  int tMaxInt;

  // parametres pour Chebychev
  int npC, npCheb, npChebev;
  double pMax, pMin;
  
  // Parametres de parallelisation
  int ompOn, numberOfThreads;

  // Parametres d'affichage
  int itSaveInLog, itSaveInData, itSaveEvthgInData; 
    
};


#endif
