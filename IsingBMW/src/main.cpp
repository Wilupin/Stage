#include "Parametres.h"
#include "Initialiseur.h"
#include "InOut.h"


// On initialise de nos objets singletons au pointeur null
Parametres *Parametres::_Parametres = NULL;
InOut      *InOut::_InOut = NULL;



int main()
{
  
  std::cout << " ------------------------------------- " << std::endl;
  std::cout << " ---------- BMW Ising Solver --------- " << std::endl;
  std::cout << " ------------------------------------- " << std::endl;
  std::cout << std::endl; 
  std::cout << "Debut de la resolution en ""temps"" : " << std::endl;

  
  // Creation des objets
  Initialiseur init;
  GammaAndDerivatives gamma;


  // Parametres d'affichage
  int itSaveInLog, itSaveInData, itSaveEvthgInData;
  itSaveInLog       = Parametres::getInstance()->get_itSaveInLog();
  itSaveInData      = Parametres::getInstance()->get_itSaveInData();
  itSaveEvthgInData = Parametres::getInstance()->get_itSaveEvthgInData();
  
  
  // Tout a ete initalisa la ou il faut on kill la classe Parametres
  Parametres::getInstance()->kill(); 
    
  
  // ========================================== //
  // Time loop
  // ========================================== //


  
  for(int itt = 0 ; itt <= gamma.get_itMax(); itt++)
    {
      if(itt==0) itt = init.Initialising(gamma);
      
      // Mise a jour de k et t
      gamma.set_kktt(itt);
      
      // Calcul des derivees et des integrales
      gamma.ChebychevDeriveesEtIntegrales();


      // Affichqge des fonctions
      if(itt%10                == 0)  std::cout << "Pas : " << itt  <<" | Temps : " << gamma.get_tt() << std::endl;
      if(itt%itSaveInLog       == 0)  gamma.SortieLogRho0();
      if(itt%itSaveInData      == 0)  gamma.SortieDataRho0();
      if(itt%itSaveEvthgInData == 0)  gamma.SortieDiversesGnuPlot(1);

      // Going one step forward in the flow equation
      gamma.FlowStepForward2D();      

    }
  
  InOut::getInstance()->CloseLog();

  return 0;
}
