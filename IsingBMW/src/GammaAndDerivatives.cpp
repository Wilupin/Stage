#include "GammaAndDerivatives.h"



// Implementation des fonctions de la classe GammaAndDerivatives



// Constructeur
GammaAndDerivatives::GammaAndDerivatives() : Der(&Cheb), GL(2)
{

  d      = 2;
  
  itMax  = Parametres::getInstance()->get_tMaxInt();
  dtt    = Parametres::getInstance()->get_dt();
  ompOn  = Parametres::getInstance()->get_ompOn();
  nrho   = Parametres::getInstance()->get_nrho();
  drho   = Parametres::getInstance()->get_drho();
  qMax   = Parametres::getInstance()->get_qMax();
  alpha  = Parametres::getInstance()->get_alpha();
  mu     = Parametres::getInstance()->get_mu();
  pMax   = Parametres::getInstance()->get_pMax();
  tMin   = Parametres::getInstance()->get_tMin();
  tMax   = Parametres::getInstance()->get_tMax();

  

  // Initialisation des variables pour le code Ising 2D
  W_rho.resize(nrho+1);
  W1_rho.resize(nrho+1);
  W2_rho.resize(nrho+1);
  V_rho.resize(nrho+1);
  
  I1_rho.resize(nrho+1);
  I2_rho.resize(nrho+1);
  I3_rho.resize(nrho+1);
  I11_rho.resize(nrho+1);

  Delta_px_py_rho.resize(Cheb.get_n());
  Delta1_px_py_rho.resize(Cheb.get_n());
  Delta2_px_py_rho.resize(Cheb.get_n());

  ChebDelta_px_py_rho.resize(Cheb.get_n());

  J3_px_py_rho.resize(Cheb.get_n());
   
  for(int ipx = 0; ipx < Cheb.get_n(); ipx++)
    {
      Delta_px_py_rho[ipx].resize(ipx+1);
      Delta1_px_py_rho[ipx].resize(ipx+1);
      Delta2_px_py_rho[ipx].resize(ipx+1);

      ChebDelta_px_py_rho[ipx].resize(ipx+1);

      J3_px_py_rho[ipx].resize(ipx+1);

      for(int ipy = 0; ipy < ipx+1; ipy++)
	{
	  Delta_px_py_rho[ipx][ipy].resize(nrho+1, 0.0);
	  Delta1_px_py_rho[ipx][ipy].resize(nrho+1, 0.0);
	  Delta2_px_py_rho[ipx][ipy].resize(nrho+1, 0.0);

	  ChebDelta_px_py_rho[ipx][ipy].resize(nrho+1, 0.0);

	  J3_px_py_rho[ipx][ipy].resize(nrho+1, 0.0);
	}
    }

 

  
  // Ecriture de la fonction de dispersion
  ptrInOut->OpenFileE0();

  double px, py; 
  
  for(int ipx = 0; ipx < Cheb.get_n(); ipx++)
    {
      px = Cheb.get_xcheb(ipx);
      
      for(int ipy = 0; ipy < Cheb.get_n(); ipy++)
	{
	  py = Cheb.get_xcheb(ipy);

	  ptrInOut->get_fileE0() << px << "\t" << py << "\t" << e0(px, py) << std::endl;
	}

    }

  ptrInOut->CloseFileE0();



  

  // Affichage de la sortie de l'initialisation

  std::cout << " | Initialisation GammaAndDerivatives complete " << std::endl;
  std::cout << " | Dimension             : " << d << std::endl;
  std::cout << " | Alpha                 : " << alpha << std::endl;
  std::cout << " | Mu                    : " << mu << std::endl; 
  std::cout << " | Temps maximal (int)   : " << itMax << std::endl;
  std::cout << " | Temps maximal (vrai)  : " << tMax << std::endl;
  std::cout << " | Temps minimal (vrai)  : " << tMin << std::endl;
  std::cout << " | Pas de temps          : " << dtt << std::endl; 
  std::cout << " | Np en rho             : " << nrho << std::endl;
  std::cout << " | Pas en rho            : " << drho << std::endl;  
  std::cout << " | Parallelisation 1/0   : " << ompOn << std::endl;
  std::cout << std::endl; 
  
}






/* -----------------------------------------------------
   Mise a jour du temps et de k ou l'on calcule
   ----------------------------------------------------- */

void GammaAndDerivatives::set_kktt(int const& n_itt)
{
  itt = n_itt;
  tt  = itt*dtt - tMin;
  kk  = exp(tt*dtt);
}

/* -------------------------------------------------
   ------------------------------------------------- *





/* -----------------------------------------------------
 Calculation of gamma derivatives and integrals 
 ------------------------------------------------------- */

void GammaAndDerivatives::ChebychevDeriveesEtIntegrales()
{
  #pragma omp parallel for if(ompOn)
  for(int ir = 0; ir <= nrho ; ir++)
    {
      
      Cheb.CoeffCheb2D(Delta_px_py_rho,
		       ChebDelta_px_py_rho, ir); // Recuperation des coeffs de Chebychev
      DerivationDim2D(ir);                       // Calcul des derivees
      IntegralsDim2D(ir);                        // Calcul des integrales
      DerivationI11Dim2D(ir);                    // Calcul de la derivee d'une integrale
    }
}





void GammaAndDerivatives::DerivationDim2D(int ir)
{
  W1_rho[ir] = Der.derivePrm(W_rho, ir, nrho);
  W2_rho[ir] = Der.deriveScd(W_rho, ir, nrho);
  
  for(int ipx = 0; ipx < Cheb.get_n(); ipx++)
    for(int ipy = 0; ipy < ipx + 1; ipy++)
    {
      // Derivatives of the functions Delta
      Delta1_px_py_rho[ipx][ipy][ir] = Der.derivePrm(Delta_px_py_rho[ipx][ipy],ir,nrho);
      Delta2_px_py_rho[ipx][ipy][ir] = Der.deriveScd(Delta_px_py_rho[ipx][ipy],ir,nrho);
    }
  
}


  
void GammaAndDerivatives::DerivationI11Dim2D(int ir)
{
  I11_rho[ir] = Der.derivePrm(I1_rho, ir, nrho);
}






/* Fonction d'integration des variables dimensionnees dans la cas Ising 2D */

void GammaAndDerivatives::IntegralsDim2D(int ir)
{

  // Variables d'integration
  double qx = 0.0, qy = 0.0; 
  
  // Integrales I
  double sumI1  = 0.0;
  double sumI2  = 0.0; 
  double sumI3  = 0.0;

  double weight = 0.0;

  // Variables pour l'integration sur p+q
  double px = 0.0, py = 0.0; 
  double pxpqx = 0.0;
  double pypqy = 0.0;
  

  // Initialisation du tableau pour J3
  for(int ipx = 0; ipx < Cheb.get_n(); ipx++)
    for(int ipy = 1; ipy < ipx+1; ipy++)
      J3_px_py_rho[ipx][ipy][ir] = 0.0;


  // Evaluation de la masse
  m2k = W_rho[ir] + 2*ir*W1_rho[ir];
  
  
  for(int j = 1; j<=GL.get_ngl(); j++)
    {
      qx = 0.5*qMax*(GL.get_x(j)+1);
      
      // On fait l'integrale sur le triangle, on se sert de la symetrie
      for (int k = 1; k <= j; k++)
	{
	  qy = 0.5*qMax*(GL.get_x(k)+1);

	  // On met a jour le propagateur 
	  UpdatePropagatorQ(qx, qy, ir);

	  // On donne les nouveaux poids
	  weight = GL.get_w(j)*GL.get_w(k); 
	  
	  if(qy < qx)
	    {
	      sumI1 = sumI1 + 2*weight*fToIntI(qx, qy, 1);
	      sumI2 = sumI2 + 2*weight*fToIntI(qx, qy, 2);
	      sumI3 = sumI3 + 2*weight*fToIntI(qx, qy, 3);
	    }
	  if(qy == qx)
	    {
	      sumI1 = sumI1 + weight*fToIntI(qx, qy, 1);
	      sumI2 = sumI2 + weight*fToIntI(qx, qy, 2);
	      sumI3 = sumI3 + weight*fToIntI(qx, qy, 3);
	    }


	  // Calcul de l'integrale J3
	  for(int ipx = 0; ipx < Cheb.get_n(); ipx ++)
	    {
	      px = Cheb.get_xcheb(ipx);
	      
	      // std::cout << j << " " << k << " " << ipx << " Bonjour" << std::endl;
	      
	      for(int ipy = 0; ipy < ipx+1; ipy++)
		{
		  py = Cheb.get_xcheb(ipy);

		  // std::cout << ipy << " Bonjour" << std::endl; 

		  Recenter(px, py, qx, qy, pxpqx, pypqy);
		  UpdatePropagatorPQ(pxpqx, pypqy, ir);

		  // std::cout << ipy << " Bonjour" << std::endl;
		  
		  if(qy < qx)
		    J3_px_py_rho[ipx][ipy][ir] += 2*weight*fToIntJ(qx, qy, 3) * GL.get_norme();
		  if(qy == qx)
		    J3_px_py_rho[ipx][ipy][ir] +=   weight*fToIntJ(qx, qy, 3) * GL.get_norme();

		  // std::cout << ipy << " Bonjour" << std::endl; 
		}
	    }
	}
      
    }

  // std::cout << "Bonjour" << std::endl;
  I1_rho[ir] = sumI1 * GL.get_norme();
  I2_rho[ir] = sumI2 * GL.get_norme();
  I3_rho[ir] = sumI3 * GL.get_norme();
}




void GammaAndDerivatives::UpdatePropagatorQ(double qx, double qy, int ir)
{
  // Recuperation de la valeur de Delta au bon point
  valDeltaQ = Cheb.chebev2D(ChebDelta_px_py_rho, qx, qy, ir);

  
  if(valDeltaQ != valDeltaQ)
    {
      if(omp_get_thread_num() == 0)
	{
	  std::cout << "ERROR, au temps : " << itt << " pour " << valDeltaQ <<  std::endl;
	  SortieDiversesGnuPlot(1);
	  exit(0);
	}
    }      

  
  // Ecriture du nouveau prpagateur
  propagatorQ = pow(  RegulatorDim(qx, qy) + e0(qx,qy) + valDeltaQ + m2k ,-1);
}



void GammaAndDerivatives::UpdatePropagatorPQ(double pxpqx, double pypqy, int ir)
{
  // Recuperation de la valeur de Delta au bon point
  valDeltaPQ = Cheb.chebev2D(ChebDelta_px_py_rho, pxpqx, pypqy, ir); 
  
  // Ecriture du nouveau prpagateur
  propagatorPQ = pow(  RegulatorDim(pxpqx, pypqy) + e0(pxpqx,pypqy) + valDeltaPQ  + m2k ,-1);
}



double GammaAndDerivatives::fToIntI(double qx, double qy, int i)
{
  return DerRegulatorDim(qx, qy)*pow(propagatorQ, i);
}


double GammaAndDerivatives::fToIntJ(double qx, double qy, int i)
{
  return DerRegulatorDim(qx, qy)*pow(propagatorQ, i-1)*propagatorPQ; 
}





double GammaAndDerivatives::RegulatorDim(double qx, double qy)
{
  ek = kk*kk*4*d*(d + mu)/(mu-2);
  
  return alpha*e0(qx,qy)*(1/(exp( e0(qx, qy) / ek ) -1));
}


double GammaAndDerivatives::DerRegulatorDim(double qx, double qy)
{
  ek = kk*kk*4*d*(d + mu)/(mu-2);
  
  return ( (2*alpha/(ek*kk))*pow(e0(qx,qy),2)*exp(e0(qx,qy)/ek) ) / pow(exp(e0(qx,qy)/ek) -1, 2) ;
}



double GammaAndDerivatives::e0(double qx, double qy)
{
  gammaq = (cos(qx) + cos(qy))/d;

  return (2*d*(d+mu))*(1-gammaq)/(d*gammaq + mu);
}


// Fonction qui utilise la periodicité pour calculer les integrales
// sur p+q la ou p+q sort du carre/triangle d'integration
void GammaAndDerivatives::Recenter(double px, double py, double qx, double qy,
				   double & npxpqx, double & npypqy)
{

  // Initialisation des coordonnes
  npxpqx = px+qx;
  npypqy = py+qy;
  
  if(px + qx > pMax)
    {
      // On revient dans la bonne zone
      npxpqx = 2*pMax - px - qx;  
    }
  if(py + qy > pMax)
    {
      npypqy = 2*pMax -py - qy;
    }
}





/* -------------------------------------------------
   ------------------------------------------------- */











/* -------------------------------------------------
   Avancee du flot pour le cas Ising 2D
   Equations de flot
   ------------------------------------------------- */

// Function called in main to make one setp forward in time
void GammaAndDerivatives::FlowStepForward2D()
{
    
#pragma omp parallel for if(ompOn)
  for(int ir = 0; ir <= nrho ; ir++)
    {
      W_rho[ir] += dtt*FlowOfWDim(ir);
      V_rho[ir] += dtt*FlowOfVDim(ir);
	
      for(int ipx = 0; ipx < Cheb.get_n(); ipx++)
	for(int ipy = 0; ipy < ipx + 1; ipy++)
	  Delta_px_py_rho[ipx][ipy][ir] += dtt*FlowOfDeltaDim(ipx, ipy, ir);
    }

}



double GammaAndDerivatives::FlowOfDeltaDim(int ipx, int ipy, int ir)
{

  uk = 2*ir*W2_rho[ir]/drho + 3*W1_rho[ir]/drho;
    
  return 2*ir*drho*J3_px_py_rho[ipx][ipy][ir] *
    pow( uk + Delta1_px_py_rho[ipx][ipy][ir]/drho, 2) - 2*ir*drho*I3_rho[ir]*pow(uk,2) -
    0.5*I2_rho[ir]*(Delta1_px_py_rho[ipx][ipy][ir] + 2*ir* Delta2_px_py_rho[ipx][ipy][ir])/drho;
}



double GammaAndDerivatives::FlowOfWDim(int ir)
{
  return 0.5*I11_rho[ir]; 
}


double GammaAndDerivatives::FlowOfVDim(int ir)
{
  return 0.5*I1_rho[ir];
}


/* -------------------------------------------------
   ------------------------------------------------- */














/* ------------------------------------------------
   Fonctions s'occupant des sorties 
   ------------------------------------------------ */


void GammaAndDerivatives::set_irho0()
{
  
  double wi=W_rho[0],wip1=W_rho[1];
  int ir = 0;
    
  while(!(wi<=0&&wip1>0)&&ir<nrho)
    {
      ir++;
      wi=W_rho[ir];
      wip1=W_rho[ir+1];
    }
    
  irho0 = ir;
  
}




void GammaAndDerivatives::SortieDataRho0()
{
  // Mise a jour de irho0
  set_irho0();

  // Ecriture de irho0
  ptrInOut->OpenFileRho0();
  ptrInOut->get_fileRho0() << tt*dtt << "\t" << irho0*drho  << std::endl; 
  ptrInOut->CloseFileRho0();
}



void GammaAndDerivatives::SortieLogRho0()
{
  // Mise a jour de irho0
  set_irho0();

  // Ecriture de irho0
  ptrInOut->get_log() << tt*dtt << "\t" << irho0*drho  << std::endl; 
}



void GammaAndDerivatives::SortieDiversesGnuPlot(int nrun)
{

  // Ouverture des fichiers de sortie
  ptrInOut->OpenSortiesGnuplot();

  for(int ir = 0; ir <= nrho; ir++)
    {
      ptrInOut->get_fileW()  << nrun << "\t" << itt << "\t"
			     << ir*drho << "\t" << W_rho[ir] << std::endl;
      ptrInOut->get_fileW1() << nrun << "\t" << itt << "\t"
			     << ir*drho << "\t" << W1_rho[ir] << std::endl;
      ptrInOut->get_fileW2() << nrun << "\t" << itt << "\t"
			     << ir*drho << "\t" << W2_rho[ir] << std::endl;
      ptrInOut->get_fileV()  << nrun << "\t" << itt << "\t" << ir*drho << "\t"
			     << V_rho[ir] << std::endl;

      ptrInOut->get_fileI1()  << nrun << "\t" << itt << "\t"
			      << ir*drho << "\t" << I1_rho[ir] << std::endl;
      ptrInOut->get_fileI2()  << nrun << "\t" << itt << "\t"
			      << ir*drho << "\t" << I2_rho[ir] << std::endl;
      ptrInOut->get_fileI3()  << nrun << "\t" << itt << "\t"
			      << ir*drho << "\t" << I3_rho[ir] << std::endl;
      ptrInOut->get_fileI11() << nrun << "\t" << itt << "\t"
			      << ir*drho << "\t" << I11_rho[ir] << std::endl;

      for(int ipx = 0; ipx < Cheb.get_n(); ipx++)
	{
	  
	  for(int ipy = 0; ipy < Cheb.get_n(); ipy++)
	    {
	      if(ipy <= ipx)
		{
		  ptrInOut->get_fileDelta() << nrun << "\t" << itt << "\t" << ir*drho << "\t"
				       << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				       << "\t" << Delta1_px_py_rho[ipx][ipy][ir] << std::endl;
		  ptrInOut->get_fileDelta1() << nrun << "\t" << itt << "\t"<<  ir*drho << "\t"
				       << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				       << "\t" << Delta1_px_py_rho[ipx][ipy][ir] << std::endl;
		  ptrInOut->get_fileDelta2() << nrun << "\t" << itt << "\t" << ir*drho << "\t"
				       << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				       << "\t" << Delta2_px_py_rho[ipx][ipy][ir] << std::endl;
		  ptrInOut->get_fileJ3() << nrun << "\t" << itt << "\t" <<  ir*drho << "\t"
				    << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				    << "\t" << J3_px_py_rho[ipx][ipy][ir] << std::endl;
		}
	      else
		{
		  ptrInOut->get_fileDelta() << nrun << "\t" << itt << "\t" << ir*drho << "\t"
				       << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				       << "\t" << Delta1_px_py_rho[ipy][ipx][ir] << std::endl;
		  ptrInOut->get_fileDelta1() << nrun << "\t" << itt << "\t" << ir*drho << "\t"
				       << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				       << "\t" << Delta1_px_py_rho[ipy][ipx][ir] << std::endl;
		  ptrInOut->get_fileDelta2() << nrun << "\t" << itt << "\t" << ir*drho << "\t"
				       << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				       << "\t" << Delta2_px_py_rho[ipy][ipx][ir] << std::endl;
		  ptrInOut->get_fileJ3() << nrun << "\t" << itt << "\t" << ir*drho << "\t"
				    << Cheb.get_xcheb(ipx) << "\t" << Cheb.get_xcheb(ipy)
				    << "\t" << J3_px_py_rho[ipy][ipx][ir] << std::endl;
		}
	    
	    }

	  ptrInOut->get_fileDelta() << std::endl;
	  ptrInOut->get_fileDelta1() << std::endl;
	  ptrInOut->get_fileDelta2() << std::endl;
	  ptrInOut->get_fileJ3() << std::endl;
	  
	}

    }
  

  ptrInOut->get_fileW() << std::endl;
  ptrInOut->get_fileW1() << std::endl;
  ptrInOut->get_fileW2() << std::endl;
  ptrInOut->get_fileV() << std::endl;

  ptrInOut->get_fileI1() << std::endl; 
  ptrInOut->get_fileI2() << std::endl;
  ptrInOut->get_fileI3() << std::endl;
  ptrInOut->get_fileI11() << std::endl;

  ptrInOut->get_fileDelta() << std::endl;
  ptrInOut->get_fileDelta1() << std::endl;
  ptrInOut->get_fileDelta2() << std::endl;
  ptrInOut->get_fileJ3() << std::endl;
     

  ptrInOut->CloseSortiesGnuplot(); 

}
