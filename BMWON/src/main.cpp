#include "initialisationEtDichotomie.h"
#include "function.h"
#include "deriveesEtEquations.h"

#define tsave 10000
#define tout 100
#define tdisplay 1000

using namespace std;

int main()
{
  GammaAndDerivatives gamma;
  GaussLegendre gausslegvar = {ngausslegendre};
  DichoParameter dichoprm;

  // Call the intialisation function
  ReadAndInitialise(gamma,gausslegvar,dichoprm);

    
  // Output and precision displayed
  ofstream log("log",ios::app);
  log.precision(16);

    
  // Time loop
  for(int ttt = 0 ; ttt <= gamma.tmax; ttt++)
    {
      if(ttt==0) ttt=initialising(gamma,dichoprm);
      gamma.kk=exp(ttt*dtt);

      // Calcul des derivees et des integrales qui
      // sont utilises a cette etape du calcul
      DeriveesEtIntegrales(gamma,gausslegvar);

      if(ttt%tdisplay==0) log<<ttt*dtt<<"\t"<<gamma.eta<<"\t"<<gamma.WWW_rho[0]<<"\t"<<(gamma.irho0=set_irho_0(gamma.WWW_rho))<<endl;
      if(ttt%tout==0) OutEtaW0(ttt,gamma.eta,gamma.WWW_rho[0],gamma.irho0);
      if(ttt%tsave==0) sortieDiversesGnuPlot(gamma,ttt,dichoprm.nrun);

      // Going on step forward in the flow equation
      flowstepfoward(gamma,log);

      if(dichoprm.OnOff&&(gamma.WWW_rho[0]>dichoprm.rht||gamma.WWW_rho[0]<dichoprm.rlt)) ttt = dichotomie(gamma.WWW_rho[0],dichoprm,ttt,log);

    }

  log<<"Paramètre dichotomie fin :\n rmin = " << dichoprm.rmin<<",\n rmax = " << dichoprm.rmax<<",\n u0 = " << dichoprm.u0<<",\n tht = " << dichoprm.tht<<",\n tlt = " << dichoprm.tlt<<",\n nrun = " << dichoprm.nrun<<endl;
  log.close();

  return 0;
}