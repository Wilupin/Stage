#include "GammaAndDerivatives.h"



// Implementation des fonctions de la classe GammaAndDerivatives



// Constructeur
GammaAndDerivatives::GammaAndDerivatives(GaussLegendre* m_gausslegvar)
{

  gausslegvar = m_gausslegvar;

  
  YYYa_p_rho.resize(nppp);
  YYYa1_p_rho.resize(nppp);
  YYYa2_p_rho.resize(nppp);
  YYYb_p_rho.resize(nppp);
  YYYb1_p_rho.resize(nppp);
  YYYb2_p_rho.resize(nppp);
  
  JJJ3LL_p_rho.resize(nppp);
  JJJ3TT_p_rho.resize(nppp);
  JJJ3LT_p_rho.resize(nppp);
  JJJ3TL_p_rho.resize(nppp);
  
  chebcoef_YYYa1_rho_p.resize(nrho+1);
  chebcoef_YYYa_rho_p.resize(nrho+1);
  chebcoef_YYYap1_rho_p.resize(nrho+1);
  chebcoef_YYYb1_rho_p.resize(nrho+1);
  chebcoef_YYYb_rho_p.resize(nrho+1);
  chebcoef_YYYbp1_rho_p.resize(nrho+1);

  III2LL_rho.resize(nrho+1,0.0);
  III2TT_rho.resize(nrho+1,0.0);
  III3LL_rho.resize(nrho+1,0.0);
  III3TT_rho.resize(nrho+1,0.0);
  III3LT_rho.resize(nrho+1,0.0);
  III3TL_rho.resize(nrho+1,0.0);
  IIIA_rho.resize(nrho+1,0.0);
  WWWI1_rho.resize(nrho+1,0);
  
  WWW_rho.resize(nrho+1,0.0);
  WWW1_rho.resize(nrho+1,0);
  WWW2_rho.resize(nrho+1,0);


  
  for(int ipp = 0; ipp < nppp; ipp++)
    {
      YYYa_p_rho[ipp].resize(nrho+1,0.0);
      YYYa1_p_rho[ipp].resize(nrho+1,0.0);
      YYYa2_p_rho[ipp].resize(nrho+1,0.0);
      YYYb_p_rho[ipp].resize(nrho+1,0.0);
      YYYb1_p_rho[ipp].resize(nrho+1,0.0);
      YYYb2_p_rho[ipp].resize(nrho+1,0.0);
      JJJ3LL_p_rho[ipp].resize(nrho+1,0.0);
      JJJ3TT_p_rho[ipp].resize(nrho+1,0.0);
      JJJ3LT_p_rho[ipp].resize(nrho+1,0.0);
      JJJ3TL_p_rho[ipp].resize(nrho+1,0.0);
    }

  for(int iii = 0; iii <= nrho; iii++)
    {
      chebcoef_YYYa1_rho_p[iii].resize(npcheb,0.0);
      chebcoef_YYYa_rho_p[iii].resize(npcheb,0.0);
      chebcoef_YYYap1_rho_p[iii].resize(npcheb,0.0);
      chebcoef_YYYb1_rho_p[iii].resize(npcheb,0.0);
      chebcoef_YYYb_rho_p[iii].resize(npcheb,0.0);
      chebcoef_YYYbp1_rho_p[iii].resize(npcheb,0.0);
    }
}




// Function called in main to make one setp forward in time
void GammaAndDerivatives::flowstepfoward(std::ofstream & log)
{
  eta = etastepforward();
    
#pragma omp parallel for if(ompOn)
  for(int iii = 0; iii <= nrho ; iii++)
    {
      WWW_rho[iii] += dtt*flowOfW(iii);
	
      for(int ipp = 0; ipp < nppp; ipp++)
        {
	  YYYa_p_rho[ipp][iii] += dtt*flowOfYa(ipp,iii);
	  YYYb_p_rho[ipp][iii] += dtt*flowOfYb(ipp,iii);
        }
    }
    
  Zk -= eta*Zk*dtt;
  YYYa_p_rho[pPrescription][rhoPrescription]=0.0;
}




double GammaAndDerivatives::flowOfYa(int ipp, int iii)
{
  
  return eta*(1+YYYa_p_rho[ipp][iii]) +
    XXX->get_xcheb(ipp)*XXX->chebev(chebcoef_YYYap1_rho_p[iii],XXX->get_xcheb(ipp)) +
    (d-2+eta)*iii*YYYa1_p_rho[ipp][iii] +
    2*iii*pow(XXX->get_xcheb(ipp),-2)*( pow(pow(XXX->get_xcheb(ipp),2)*YYYa1_p_rho[ipp][iii] +
					    WWW1_rho[iii],2)*JJJ3LT_p_rho[ipp][iii] +
					pow(pow(XXX->get_xcheb(ipp),2)*YYYb_p_rho[ipp][iii]*drho + WWW1_rho[iii],2)*JJJ3TL_p_rho[ipp][iii] -
					pow(WWW1_rho[iii],2)*(III3LT_rho[iii]+III3TL_rho[iii]) )/drho -
    0.5*III2LL_rho[iii]*(YYYa1_p_rho[ipp][iii]+2*iii*YYYa2_p_rho[ipp][iii])/drho -
    0.5*III2TT_rho[iii]*((n-1)*YYYa1_p_rho[ipp][iii]/drho + 2*YYYb_p_rho[ipp][iii]);
}





double GammaAndDerivatives::flowOfYb(int ipp, int iii)
{
    
  return (2*eta-2+d)*YYYb_p_rho[ipp][iii] +
    XXX->get_xcheb(ipp)*XXX->chebev(chebcoef_YYYbp1_rho_p[iii],XXX->get_xcheb(ipp)) +
    (d-2+eta)*iii*YYYb1_p_rho[ipp][iii] +
    ( (n-1)*(JJJ3TT_p_rho[ipp][iii]*pow(pow(XXX->get_xcheb(ipp),2)*YYYb_p_rho[ipp][iii] + WWW1_rho[iii]/drho,2) -
	     III3TT_rho[iii]*pow(WWW1_rho[iii]/drho,2)) -
      ( JJJ3LT_p_rho[ipp][iii]*pow(pow(XXX->get_xcheb(ipp),2)*YYYa1_p_rho[ipp][iii]/drho + WWW1_rho[iii]/drho,2) +
	JJJ3TL_p_rho[ipp][iii]*pow(pow(XXX->get_xcheb(ipp),2)*YYYb_p_rho[ipp][iii] + WWW1_rho[iii]/drho,2) -
	(III3LT_rho[iii]+III3TL_rho[iii])*pow(WWW1_rho[iii]/drho,2) ) +
      JJJ3LL_p_rho[ipp][iii]*pow(pow(XXX->get_xcheb(ipp),2)*(YYYa1_p_rho[ipp][iii]/drho + 2*YYYb_p_rho[ipp][iii] +
							     2*iii*YYYb1_p_rho[ipp][iii]) + 3*WWW1_rho[iii]/drho +
				 2*iii*WWW2_rho[iii]/drho,2) - pow(3*WWW1_rho[iii]/drho +
								   2*iii*WWW2_rho[iii]/drho,2)*III3LL_rho[iii] )*pow(XXX->get_xcheb(ipp),-2) -
    0.5*III2LL_rho[iii]*(5*YYYb1_p_rho[ipp][iii]+2*iii*YYYb2_p_rho[ipp][iii])/drho -
    0.5*III2TT_rho[iii]*(n-1)*YYYb1_p_rho[ipp][iii]/drho + IIIA_rho[iii]*YYYb_p_rho[ipp][iii];
}




double GammaAndDerivatives::flowOfW(int iii)
{
  return (eta-2)*WWW_rho[iii] + (d - 2 + eta)*iii*WWW1_rho[iii] + 0.5*WWWI1_rho[iii];
}



double GammaAndDerivatives::etastepforward()
{
  int iii = rhoPrescription, ipp = pPrescription;
  
  return -(XXX->get_xcheb(ipp)*XXX->chebev(chebcoef_YYYap1_rho_p[iii],XXX->get_xcheb(ipp)) +
	   (d-2)*iii*YYYa1_p_rho[ipp][iii] +
	   2*iii*pow(XXX->get_xcheb(ipp),-2)*( pow(pow(XXX->get_xcheb(ipp),2)*YYYa1_p_rho[ipp][iii] +
						   WWW1_rho[iii],2)*JJJ3LT_p_rho[ipp][iii] +
					       pow(pow(XXX->get_xcheb(ipp),2)*YYYb_p_rho[ipp][iii]*drho +
						   WWW1_rho[iii],2)*JJJ3TL_p_rho[ipp][iii] -
					       pow(WWW1_rho[iii],2)*(III3LT_rho[iii]+III3TL_rho[iii]) )/drho -
	   0.5*III2LL_rho[iii]*(YYYa1_p_rho[ipp][iii]+2*iii*YYYa2_p_rho[ipp][iii])/drho -
	   0.5*III2TT_rho[iii]*((n-1)*YYYa1_p_rho[ipp][iii]/drho + 2*YYYb_p_rho[ipp][iii]))/(1+iii*YYYa1_p_rho[ipp][iii]);
}











// Calculation of Chebytchev coeff for the gamma functions
void GammaAndDerivatives::CoeffChebgamma()
{
    //On actualise les coefficients chebyshev
    #pragma omp parallel for if(ompOn)
    for(int iii = 0 ; iii <= nrho; iii++)
      {
	XXX->CoeffCheb2(YYYa_p_rho, chebcoef_YYYa_rho_p, iii);
        XXX->CoeffCheb2(YYYb_p_rho, chebcoef_YYYb_rho_p, iii);
      }
}




// Calculation of derivatives in "rho = rho[iii]" and with the spectral method in p
void GammaAndDerivatives::Derivation(int iii)
{
  
  // Calculation of Chebytchev coefficients for YYYa_rho_p and YYYb_rho_p
  XXX->multiCoeffCheb2_2(YYYa_p_rho, chebcoef_YYYa_rho_p,YYYb_p_rho, chebcoef_YYYb_rho_p,iii);

  // Derivative of YYYa_rho_p and YYYb_rho_p with respect to p (using Chebytchev polynoms)
  Der->multichder_2(chebcoef_YYYa_rho_p[iii],chebcoef_YYYap1_rho_p[iii],
		    chebcoef_YYYb_rho_p[iii],chebcoef_YYYbp1_rho_p[iii]);


  // Derivative of WWW1_rho and WWW2_rho (which does not depend on p) with respect to rho
  WWW1_rho[iii] = Der->derivePrm(WWW_rho,iii,nrho);
  WWW2_rho[iii] = Der->deriveScd(WWW_rho,iii,nrho);
  
  // Derivative with respect to rho for each given p of the other functions 
  for(int ipp = 0; ipp < nppp; ipp++)
    {
      // Derivatives of the integrals
      JJJ3LL_p_rho[ipp][iii] = 0.0;
      JJJ3TT_p_rho[ipp][iii] = 0.0;
      JJJ3LT_p_rho[ipp][iii] = 0.0;
      JJJ3TL_p_rho[ipp][iii] = 0.0;

      // Derivatives of the functions
      YYYa1_p_rho[ipp][iii] = Der->derivePrm(YYYa_p_rho[ipp],iii,nrho);
      YYYa2_p_rho[ipp][iii] = Der->deriveScd(YYYa_p_rho[ipp],iii,nrho);
      YYYb1_p_rho[ipp][iii] = Der->derivePrm(YYYb_p_rho[ipp],iii,nrho);
      YYYb2_p_rho[ipp][iii] = Der->deriveScd(YYYb_p_rho[ipp],iii,nrho);
    }
  
  XXX->multiCoeffCheb2_2(YYYa1_p_rho, chebcoef_YYYa1_rho_p,YYYb1_p_rho, chebcoef_YYYb1_rho_p,iii);

}





// Integrals calculation
void GammaAndDerivatives::Integrals(int iii)
{
  
  double bpaq1 = 0, bmaq1 = qInt, bpaq2 = 0.5*qInt, bmaq2 = 0.5*qInt;
  double q1, q2, p, qsqr, facReg, propgl, propgt, auxl, auxt, yya, yyb, yyb1, yya1;

  double mk2    = WWW_rho[iii] + 2*iii*WWW1_rho[iii]; // WARNING !!
  double www    = WWW_rho[iii];
  double www1   = WWW1_rho[iii];
    
  double int2ll = 0.0, int2tt = 0.0, int3ll = 0.0;
  double int3tt = 0.0, int3lt = 0.0, int3tl = 0.0;
  double inta   = 0.0, int1p = 0.0;
    
  double yamax, ybmax,ya_p,yb_p;

  // Value of Ya and Yb en rho(iii) et pmax
  XXX->multichebev_2(chebcoef_YYYa_rho_p[iii],chebcoef_YYYb_rho_p[iii],yamax,ybmax,p_max);


  // Main integration loop
  for(int jjj=1;jjj<=gausslegvar->get_ngl();jjj++)
    {
      // Variable de la premiere integrale
      q2 = (bpaq2+bmaq2*gausslegvar->get_x(jjj));
      
      for(int kkk=1;kkk<=gausslegvar->get_ngl();kkk++)
	{
	  // Variable de la seconde integrale
	  q1 = (bpaq1+bmaq1*gausslegvar->get_x(kkk));
	  qsqr = pow(q1,2)+pow(q2,2);

	  facReg = gausslegvar->get_w(jjj)*gausslegvar->get_w(kkk)*pow(q2,gausslegvar->get_d()-2)*
	    sss(qsqr,gausslegvar->get_alpha(),eta);

	  XXX->multichebev_4(chebcoef_YYYa_rho_p[iii],chebcoef_YYYb_rho_p[iii],
				    chebcoef_YYYa1_rho_p[iii],chebcoef_YYYb1_rho_p[iii],
				    yya,yyb,yya1,yyb1,sqrt(qsqr));

	  // Propagators (Longitudinal and transversal)
	  propgl = pow( regulator(qsqr,gausslegvar->get_alpha()) +
			qsqr*( 1 + yya + 2*iii*drho*yyb ) + mk2,-1 ) ;
	  propgt = pow( regulator(qsqr,gausslegvar->get_alpha()) +
			qsqr*( 1 + yya ) + www,-1 ) ;

	  auxl = facReg*pow(propgl,2);
	  auxt = facReg*pow(propgt,2);

	  int1p += -(n-1)*auxt*(qsqr*yya1 + www1)/drho -
	    auxl*(qsqr*( yya1/drho + 2*yyb + 2*iii*yyb1 ) +
		  (3*www1 + 2*iii*WWW2_rho[iii])/drho);
	  
	  inta += facReg*propgl*propgt*(propgl+propgt)*(qsqr*yyb + www1/drho);
	  int2ll += auxl;
	  int2tt += auxt;
	  int3ll += auxl * propgl;
	  int3tt += auxt * propgt;
	  int3lt += auxl * propgt;
	  int3tl += auxt * propgl;
	  
	  
	  // Loop to compute the integrals J depending on p
	  for(int ipp=0;ipp<nppp;ipp++)
	    {
	      p = XXX->get_xcheb(ipp);
	      qsqr = pow(q1,2) + pow(q2,2) + pow(p,2) + 2*p*q1;
	      p = sqrt(qsqr);
	      
	      if(p<p_max)
                {
		  XXX->multichebev_2(chebcoef_YYYa_rho_p[iii],chebcoef_YYYb_rho_p[iii],ya_p,yb_p,p);
		  propgl = pow(regulator(qsqr,gausslegvar->get_alpha()) +
			       qsqr*(1 + ya_p + 2*iii*drho*yb_p)+mk2,-1);
		  propgt = pow( regulator(qsqr,gausslegvar->get_alpha()) + qsqr*(1 + ya_p) + www,-1);
                }
	      else // if p >= p_max
		{
		  propgl = pow( regulator(qsqr,gausslegvar->get_alpha()) +
				qsqr*( 1 + yamax + 2*iii*drho*ybmax ) + mk2,-1);
		  propgt = pow( regulator(qsqr,gausslegvar->get_alpha()) + qsqr*(1 + yamax) + www,-1);
		}
	      JJJ3LL_p_rho[ipp][iii] += gausslegvar->get_norme() * auxl * propgl;
	      JJJ3TT_p_rho[ipp][iii] += gausslegvar->get_norme() * auxt * propgt;
	      JJJ3LT_p_rho[ipp][iii] += gausslegvar->get_norme() * auxl * propgt;
	      JJJ3TL_p_rho[ipp][iii] += gausslegvar->get_norme() * auxt * propgl;
            }
        }
    }

  // Finalization of the calculation on integrals I
  WWWI1_rho[iii] = int1p * gausslegvar->get_norme();
  IIIA_rho[iii] = inta * gausslegvar->get_norme();
  III2LL_rho[iii] = int2ll * gausslegvar->get_norme();
  III2TT_rho[iii] = int2tt * gausslegvar->get_norme();
  III3LL_rho[iii] = int3ll * gausslegvar->get_norme();
  III3TT_rho[iii] = int3tt * gausslegvar->get_norme();
  III3LT_rho[iii] = int3lt * gausslegvar->get_norme();
  III3TL_rho[iii] = int3tl * gausslegvar->get_norme();
}





// Calculation of gamma derivatives and integrals using previous functions
void GammaAndDerivatives::DeriveesEtIntegrales()
{
  #pragma omp parallel for if(ompOn)
  for(int iii = 0; iii <= nrho ; iii++) Derivation(iii), Integrals(iii);
}



double GammaAndDerivatives::regulator(double y, double alpha)
{
    return alpha*y/(exp(y)-1);
}

double GammaAndDerivatives::sss(double y, double alpha, double eta)
{
    return alpha*y*(-eta+2*y*exp(y)/(exp(y)-1))/(exp(y)-1);
}














void GammaAndDerivatives::sortieDiversesGnuPlot(int ttt, int nrun)
{
  
  std::ofstream outfffa("../data/fffa",std::ios::app);
  std::ofstream outfffa1("../data/fffa1",std::ios::app);
  std::ofstream outfffa2("../data/fffa2",std::ios::app);
  std::ofstream outfffap1("../data/fffap1",std::ios::app);
  std::ofstream outfffb("../data/fffb",std::ios::app);
  std::ofstream outfffb1("../data/fffb1",std::ios::app);
  std::ofstream outfffb2("../data/fffb2",std::ios::app);
  std::ofstream outfffbp1("../data/fffbp1",std::ios::app);
  std::ofstream outWWW("../data/W",std::ios::app);
  std::ofstream outWWW1("../data/W1",std::ios::app);
  std::ofstream outWWW2("../data/W2",std::ios::app);
  std::ofstream outdI1("../data/dI1",std::ios::app);
  std::ofstream outI2ll("../data/I2ll",std::ios::app);
  std::ofstream outI2tt("../data/I2tt",std::ios::app);
  std::ofstream outI3ll("../data/I3ll",std::ios::app);
  std::ofstream outI3tt("../data/I3tt",std::ios::app);
  std::ofstream outI3lt("../data/I3lt",std::ios::app);
  std::ofstream outI3tl("../data/I3tl",std::ios::app);
  std::ofstream outJ3ll("../data/J3ll",std::ios::app);
  std::ofstream outJ3tt("../data/J3tt",std::ios::app);
  std::ofstream outJ3lt("../data/J3lt",std::ios::app);
  std::ofstream outJ3tl("../data/J3tl",std::ios::app);
  std::ofstream outIa("../data/Ia");
  
  outfffa.precision(16),outfffa1.precision(16),outfffa2.precision(16);
  outfffap1.precision(16),outfffb.precision(16),outfffb1.precision(16);
  outfffb2.precision(16),outfffbp1.precision(16),outWWW.precision(16);
  outWWW1.precision(16),outWWW2.precision(16),outdI1.precision(16);
  outI2ll.precision(16),outI2tt.precision(16),outI3ll.precision(16);
  outI3tt.precision(16),outI3lt.precision(16),outI3tl.precision(16);
  outJ3ll.precision(16),outJ3tt.precision(16),outJ3lt.precision(16);
  outJ3tl.precision(16),outIa.precision(16);

  for(int iii=0; iii <= nrho; iii++)
    {
      outdI1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<WWWI1_rho[iii]<<std::endl;
      outIa<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<IIIA_rho[iii]<<std::endl;
      outI2ll<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<III2LL_rho[iii]<<std::endl;
      outI2tt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<III2TT_rho[iii]<<std::endl;
      outI3ll<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<III3LL_rho[iii]<<std::endl;
      outI3tt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<III3TT_rho[iii]<<std::endl;
      outI3lt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<III3LT_rho[iii]<<std::endl;
      outI3tl<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<III3TL_rho[iii]<<std::endl;
      outWWW <<nrun<<"\t"<<ttt<<"\t"<< iii*drho<<"\t"<<WWW_rho[iii]<<std::endl;
      outWWW1 <<nrun<<"\t"<<ttt<<"\t"<< iii*drho<<"\t"<<WWW1_rho[iii]<<std::endl;
      outWWW2 <<nrun<<"\t"<<ttt<<"\t"<< iii*drho<<"\t"<<WWW2_rho[iii]<<std::endl;

      for(int ipp=0 ; ipp < nppp ; ipp++)
        {
	  
	  outfffa<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<YYYa_p_rho[ipp][iii]<<std::endl;
	  outfffa1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<YYYa1_p_rho[ipp][iii]<<std::endl;
	  outfffa2<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<YYYa2_p_rho[ipp][iii]<<std::endl;
														   outfffap1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<XXX->chebev(chebcoef_YYYap1_rho_p[iii],XXX->get_xcheb(ipp))<<std::endl;
																			  outfffb<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<YYYb_p_rho[ipp][iii]<<std::endl;
														 outfffb1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<YYYb1_p_rho[ipp][iii]<<std::endl;
																											  outfffb2<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<YYYb2_p_rho[ipp][iii]<<std::endl;
														   outfffbp1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<XXX->chebev(chebcoef_YYYbp1_rho_p[iii], XXX->get_xcheb(ipp))<<std::endl;
	  outJ3ll<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<JJJ3LL_p_rho[ipp][iii]<<std::endl;
	  outJ3tt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<JJJ3TT_p_rho[ipp][iii]<<std::endl;
	  outJ3lt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<JJJ3LT_p_rho[ipp][iii]<<std::endl;
	  outJ3tl<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<XXX->get_xcheb(ipp)<<"\t"<<JJJ3TL_p_rho[ipp][iii]<<std::endl;
        }
      outfffa<<std::endl;
      outfffa1<<std::endl;
      outfffa2<<std::endl;
      outfffap1<<std::endl;
      outfffb<<std::endl;
      outfffb1<<std::endl;
      outfffb2<<std::endl;
      outfffbp1<<std::endl;
      outJ3ll<<std::endl;
      outJ3tt<<std::endl;
      outJ3lt<<std::endl;
      outJ3tl<<std::endl;
    }

  outfffa.close(),outfffa1.close(),outfffa2.close();
  outfffap1.close(),outfffb.close(),outfffb1.close();
  outfffb2.close(),outfffbp1.close(),outWWW.close();
  outWWW1.close(),outWWW2.close(),outdI1.close();
  outI2ll.close(),outI2tt.close(),outI3ll.close();
  outI3tt.close(),outI3lt.close(),outI3tl.close();
  outIa.close(),outJ3ll.close(),outJ3tt.close();
  outJ3lt.close(),outJ3tl.close();

}


