#ifndef DEF_INITIALISEUR
#define DEF_INITIALISEUR


#include "GammaAndDerivatives.h"
#include "Parametres.h"

class Initialiseur
{
    
public:
    
  // Constructeur et destructeur
  Initialiseur();
  ~Initialiseur(){};

    
  // Fonction de classe
  int Initialising(GammaAndDerivatives& gamma);
  int Dichotomie(double w0, int ttt);

  // Pointeur vers les fichiers de sortie
  InOut* ptrInOut = InOut::getInstance();

 private:

  double beta, d, rho, mu; 
  int nCheb;
};


#endif
