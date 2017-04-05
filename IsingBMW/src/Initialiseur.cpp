#include "Initialiseur.h"




// Constructeur
Initialiseur::Initialiseur()
{

  // On initialise les parametres
  Parametres::getInstance();
  
  // Erase preivous data in data directory
  system("rm -f ../data/*");

  // Initialisation du nombre de threads
  omp_set_num_threads(Parametres::getInstance()->get_numberOfThreads());


}







int Initialiseur::Initialising(GammaAndDerivatives& gamma)
{
  /*

  // On fixe/refixe la condition initiale sur les grandeurs scalaires
  gamma.set_eta(0.0);
  gamma.set_Zk(1.0);
  gamma.set_kk(1.0);

  double u0 = 0.02;
  double r0 = 0.01;
  

  // On fixe/refixe la condition initiale sur les fonctions
  for(int iii = 0; iii <= gamma.get_nrho(); iii++)
    {
		
      gamma.set_WWW_rho(r0 + u0*iii*gamma.get_drho(), iii);
	  
      for(int ipp = 0; ipp < gamma.get_nppp(); ipp++)
	{
	  gamma.set_YYYa_p_rho(0.0, ipp, iii);
	  gamma.set_YYYb_p_rho(0.0, ipp, iii);
	}
    }
  */
      
  return 0;
    
}








