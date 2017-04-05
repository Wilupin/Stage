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

    
};


#endif
