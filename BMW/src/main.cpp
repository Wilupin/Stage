#include "initialisationEtDichotomie.h"

#define tsave 10000
#define tout 100
#define tdisplay 1000


// On initialise notre unique objet Chebychev a null
Chebychev *Chebychev::_singletonChebychev = NULL;
Derivateur *Derivateur::_singletonDerivateur = NULL;


int main()
{
  std::cout << " " << std::endl;
  std::cout << "----------- BMW Solver ----------" << std::endl;
  std::cout << " " << std::endl;

  GaussLegendre gausslegvar;
  GammaAndDerivatives gamma(&gausslegvar);
  DichoParameter dichoprm;

  // Call the intialisation function
  ReadAndInitialise(gamma,gausslegvar,dichoprm);
    
  // Output and precision displayed
  ofstream log("log",ios::app);
  log.precision(16);


  std::cout << "Debut de la resolution en ""temps"" :" << std::endl;



  // Time loop
  for(int ttt = 0 ; ttt <= gamma.get_tmax(); ttt++)
    {
      if(ttt%10==0) std::cout << "Pas : " << ttt  <<" | Temps : " << ttt*dtt << std::endl;
      
      if(ttt==0) ttt=initialising(gamma,dichoprm);
      gamma.set_kk(exp(ttt*dtt));
      
      // Calcul des derivees et des integrales qui
      // sont utilises a cette etape du calcul
      gamma.DeriveesEtIntegrales();

      if(ttt%tdisplay==0)
	{
	  gamma.set_irho0(set_irho_0(gamma.get_WWW_rho_vec()));
	  log<<ttt*dtt<<"\t"<<gamma.get_eta()<<"\t"<<gamma.get_WWW_rho(0)<<"\t"<<gamma.get_irho0()<<endl;
	}
      if(ttt%tout==0) OutEtaW0(ttt,gamma.get_eta(),gamma.get_WWW_rho(0),gamma.get_irho0());
      if(ttt%tsave==0) gamma.sortieDiversesGnuPlot(ttt,dichoprm.nrun);

      // Going on step forward in the flow equation
      gamma.flowstepfoward(log);

      if(dichoprm.OnOff && (gamma.get_WWW_rho(0)>dichoprm.rht||gamma.get_WWW_rho(0)<dichoprm.rlt)) ttt = dichotomie(gamma.get_WWW_rho(0),dichoprm,ttt,log);
    }

  log << "Paramètre dichotomie fin :\n rmin = " << dichoprm.rmin
      << ",\n rmax = " << dichoprm.rmax
      << ",\n u0 = " << dichoprm.u0
      << ",\n tht = " << dichoprm.tht
      << ",\n tlt = " << dichoprm.tlt
      << ",\n nrun = " << dichoprm.nrun<<endl;
  log.close();

  return 0;
}
