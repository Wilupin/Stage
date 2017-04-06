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

  // Initialisation des constantes et variables
  d = 2.0;

  // On recupere le nombre de point d'interpolation pour intialiser sur tous les points
  nCheb = Parametres::getInstance()->get_npCheb();
  mu    = Parametres::getInstance()->get_mu();
  beta  = Parametres::getInstance()->get_beta();

  std::cout << " | Initialisation Generale complete " << std::endl;
  std::cout << " | Beta                  : " << beta << std::endl;
  std::cout << " | Mu                    : " << mu  << std::endl;
  std::cout << std::endl; 
  


}







int Initialiseur::Initialising(GammaAndDerivatives& gamma)
{
  

  // On fixe/refixe la condition initiale sur les grandeurs scalaires
  gamma.set_kktt(0);

  
  // On fixe/refixe la condition initiale sur les fonctions
  for(int ir = 0; ir <= gamma.get_nrho(); ir++)
    {
      rho = ir*gamma.get_drho();
      
      // Initialisation du potentiel
      gamma.set_V_rho( (d+mu)*2*rho - log(cosh( 2*sqrt(beta)*(d+mu)*sqrt(2*rho)   )), ir );

      // Initialisation de la derivee du potentiel en evitant la singularite en 0
      if(ir == 0)
	gamma.set_W_rho( 2*(d+mu) * (1 - 2*beta*(d+mu)) , ir); 
      else
	gamma.set_W_rho( 2*(d+mu) - sqrt(2*beta/rho) * (d+mu) * tanh( 2*sqrt(beta)*(d+mu)*sqrt(2*rho)), ir);

      // Initialisation de Delta
      for(int ipx = 0; ipx < nCheb; ipx++)
	for(int ipy = 0; ipy <= ipx; ipy++)
	  gamma.set_Delta_px_py_rho(0.0, ipx, ipy, ir);
    }
  
      
  return 0;
    
}








int Initialiseur::Dichotomie(double w0, int ttt)
{

  // On saute une ligne dans les fichiers de plot
  // Ceci pour ne pas ques deux run soient relies par des points
  ptrInOut->OpenFileRho0();
  ptrInOut->get_fileRho0() << std::endl; 
  ptrInOut->CloseFileRho0();

  // Condition pour effectuer la dichotomie et stopper le run :
  /* Si jamais on se retrouve avec w0 plus grand qu'une certaine valeur
     c'est que l'on se trouve dans la phase symetrique et c'est sur.
     Si jamais on a w0 plus petit qu'une certaine valeur c'est que l'on 
     se trouve dans la phase brisee et c'est sur aussi. Il faut ajuster ces
     deux valeurs rlt et rht a la main pour ne pas que les calcul soient
     trop long, tout en faisant en sorte que l'on ne les atteignent pas dans un 
     cas ou l'on ne saurait pas vraiment dans la bonne phase */

  /* ATTENTION : Dans le cas ou ttt>tht ou ttt>tlt cela signifie que l'on a une 
     divergence de W(0) qui est plus rapide sur le run qui ce finit ici, que sur
     le run precedent. Ceci n'est pas normal car on est suppose se rapprocher 
     du point critiaque et donc rester de plus en plus longtemps dans le voisinage
     du point fixe. */

  /*
    if(sb<rlt)
    {
        if(ttt>tlt)
        {
            log << "Dichotomie à t = "<<ttt*dtt<<std::endl
	       << "Phase basse temperature rmin -> r0 "<<std::endl;
            tlt = ttt;
            rmin = r0;
        }
        else
        {
            log << "La dichotomie n'avance plus, fin du programme!"<<std::endl;
            log << "Paramètre dichotomie fin :\n rmin = " << rmin
            << ",\n rmax = " << rmax<<",\n u0 = " << u0<<",\n tht = "
            << tht<<",\n tlt = " << tlt<<std::endl;
            log.close();
            exit(0);
        }
    }
    else
    {
        if(ttt>tht)
        {
            log<<"Dichotomie à t = "<<ttt*dtt<<std::endl<<"Phase haute temperature rmax -> r0"<<std::endl;
            tht = ttt;
            rmax = r0;
        }
        else
        {
            log <<"La dichotomie n'avance plus, fin du programme!"<<std::endl;
            log <<"Paramètre dichotomie fin :\n rmin = " << rmin<<",\n rmax = "
            << rmax<<",\n u0 = " << u0 <<",\n tht = "
		<< tht<<",\n tlt = " << tlt<<std::endl;
            log.close();
            exit(0);
        }
    }

    
    /* Si tout c'est bien passe on retourne -1 pour reinitialiser le temps a cette valeur
       et relancer a 0 la boucle sur le temps. */

    return -1;
}





