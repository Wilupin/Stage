#include "GaussLegendre.h"


 
// Calculation of Gauss Legendre weight and abscissas
void GaussLegendre::gaulegWeightAbscissas(GaussLegendre& gausslegvar)
{
  
  int m,i,j;
  double z1,z,pp,p3,p2,p1; //High precision is a good idea for this routine.

  gausslegvar.x.resize(gausslegvar.ngl+2);
  gausslegvar.w.resize(gausslegvar.ngl+2);

    // The roots are symmetric in the interval, so we only have to find half of them
  m=(gausslegvar.ngl+1)/2;

  
  for(i=1;i<=m;i++)
    // Loop over the desired roots
    {
      // Starting with the above approximation the ith root,
      // we enter the main loop of refinement by Newtow's method
      z = cos(M_PI*(i-0.25)/(gausslegvar.ngl+0.5));
      do
        {
	  p1=1.0;
	  p2=0.0;

	  // Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
	  for(j=1;j<=gausslegvar.ngl;j++)
            {
	      p3=p2;
	      p2=p1;
	      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
	  //p1 is now the desired Legendre polynomial.
	  //We next compute pp, its derivative, by a standard relation
	  //involving also p2, the polynomial of one lower order.
	  pp=gausslegvar.ngl*(z*p1-p2)/(z*z-1.0);
	  z1=z;
	  //Newthon's method.
	  z=z1-p1/pp;       
        }while(fabs(z-z1)>EPS);
      
      //Scale the root to the desired interval, and put in its symmetric counterpart.
      gausslegvar.x[i]=-z;                       
      gausslegvar.x[gausslegvar.ngl+1-i]=z;
      
      //Compute the weight and its symmetric counterpart
      gausslegvar.w[i]=2./((1.0-z*z)*pp*pp);      
      gausslegvar.w[gausslegvar.ngl+1-i]=gausslegvar.w[i];
    }
  gausslegvar.norme = normalisationIntegrale(gausslegvar.d);
}



double GaussLegendre::normalisationIntegrale(double const& d)
{
  return qInt*0.5*qInt*d*tgamma(d/2.)/tgamma((d-1)/2.)/sqrt(M_PI);
}


double GaussLegendre::SetKd(double const& d)
{
    return pow(2,1-d)*pow(M_PI,-d/2)*pow(d*tgamma(d/2),-1);
}
