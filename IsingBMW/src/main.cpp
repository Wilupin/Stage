#include "Parametres.h"
#include "Initialiseur.h"
#include "InOut.h"

#define tsave 10000
#define tout 100
#define tdisplay 1000


// On initialise de nos objets singletons au pointeur null
Parametres *Parametres::_Parametres = NULL;
InOut      *InOut::_InOut = NULL;


int main()
{
  std::cout << " ------------------------------------- " << std::endl;
  std::cout << " ---------- BMW Ising Solver --------- " << std::endl;
  std::cout << " ------------------------------------- " << std::endl;
  std::cout << std::endl; 

    
    // Creation des objets
  GammaAndDerivatives gamma;
  Initialiseur init;

  // Tout a ete initalisa la ou il faut on kill la classe parametres
  Parametres::getInstance()->kill(); 
    
  // Output and precision displayed
  std::ofstream log("../output/log");
  log.precision(16);


  std::cout << "Debut de la resolution en ""temps"" : " << std::endl;

  InOut::getInstance()->get_log() << "Bonjour " << std::endl;
  
  
  // ========================================== //
  // Time loop
  // ========================================== //


  /*
  for(int ttt = 0 ; ttt <= gamma.get_tmax(); ttt++)
    {
      if(ttt%10==0)
	std::cout << "Pas : " << ttt  <<" | Temps : " << ttt*gamma.get_dt() << std::endl;
      
      if(ttt==0)
          ttt = init.Initialising(gamma);
        
      gamma.set_kk(exp(ttt*gamma.get_dt()));
      
      // Calcul des derivees et des integrales
      gamma.DeriveesEtIntegrales();

      
      if(ttt%tdisplay==0)
      {
          gamma.set_irho0();
          log<<ttt*gamma.get_dt()<<"\t"<<"\t"<<gamma.get_W_rho(0)<<"\t"<<gamma.get_irho0()<<std::endl;
      }

      
      if(ttt%tout==0)
	InOut::getInstance()->OutEtaW0(ttt,gamma.get_eta(),gamma.get_WWW_rho(0),gamma.get_irho0(), gamma.get_dt());
        
      if(ttt%tsave==0)
          gamma.sortieDiversesGnuPlot(ttt,1);

      // Going one step forward in the flow equation
      gamma.flowstepfoward(log);
      

    }


  
  log.close();

  */

  return 0;
}
