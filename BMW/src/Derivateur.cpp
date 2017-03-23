#include "Derivateur.h"


//1st Derivative
double Derivateur::der1bulk(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return (fm2/4. - 2*fm1 + 2*fp1 - fp2/4.)/3.;
}
double Derivateur::der1FstPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return -25*fm2/12. + 4*fm1 - 3*f0 + 4*fp1/3. - fp2/4.;
}
double Derivateur::der1ScdPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return (-fm2/2. - 5*fm1/3. + 3*f0 - fp1 + fp2/6.)/2.;
}
double Derivateur::der1LastPtm1(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return (fp2/2. + 5*fp1/3. - 3*f0 + fm1 - fm2/6.)/2.;
}
double Derivateur::der1LastPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return 25*fp2/12. - 4*fp1 + 3*f0 - 4*fm1/3. + fm2/4. ;
}



//2nd Derivative
double Derivateur::der2bulk(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return (-fm2/4. + 4*fm1 + 4*fp1 - fp2/4.)/3. - 5*f0/2.;
}
double Derivateur::der2FstPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return 35*fm2/12. - 26*fm1/3. + 19*f0/2. - 14*fp1/3. + 11*fp2/12.;
}
double Derivateur::der2ScdPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return (11*fm2)/12. - (5*fm1)/3. + f0/2. + fp1/3. - fp2/12.;
}
double Derivateur::der2LastPtm1(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return -fm2/12. + fm1/3. + f0/2. - (5*fp1)/3. + (11*fp2)/12.;
}
double Derivateur::der2LastPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
  return 11*fm2/12. - 14*fm1/3. + 19*f0/2. - 26*fp1/3. + 35*fp2/12.;
}




void Derivateur::chder(std::vector<double> const& c, std::vector<double> & cder)
//Given a,b,c[0..n-1], as output from routine chebft §5.8, and given n, the desired degree of approximation (length of c to be used), this routine returns the array cder[0..n-1], the Chebyshev coefficients of the derivative of the function whose coefficients are c.
{
  int j;
  double con;

  Chebychev *XXX = Chebychev::getInstance();
    
  cder[XXX->get_n()-1]=0.0; //n-1 and n-2 are special cases.
  cder[XXX->get_n()-2]=2*(XXX->get_n()-1)*c[XXX->get_n()-1];
  for (j=XXX->get_n()-3;j>=0;j--)
    cder[j]=cder[j+2]+2*(j+1)*c[j+1];// Equation (5.9.2).
  con=2.0/(XXX->get_xmax()-XXX->get_xmin());
  for (j=1;j<XXX->get_n();j++) //Normalize to the interval b-a.
    {
      cder[j] *= con;
    }
  cder[0] *= 1.0/(XXX->get_xmax()-XXX->get_xmin());
}




void Derivateur::multichder_2(std::vector<double> const& c1, std::vector<double> & cder1,
			      std::vector<double> const& c2, std::vector<double> & cder2)
//Given a,b,c[0..n-1], as output from routine chebft §5.8, and given n, the desired degree of approximation (length of c to be used),
//this routine returns the array cder[0..n-1], the Chebyshev coefficients of the derivative of the function whose coefficients are c.
{
  int j;
  double con;

  Chebychev *XXX = Chebychev::getInstance();

  int N_cheb = XXX->get_n();
    
  cder1[N_cheb-1]=0.0;
  cder2[N_cheb-1]=0.0; //n-1 and n-2 are special cases.
  cder1[N_cheb-2]=2*(N_cheb-1)*c1[N_cheb-1];
  cder2[N_cheb-2]=2*(N_cheb-1)*c2[N_cheb-1];
  for (j= N_cheb-3; j>=0; j--)
    {
      cder1[j]=cder1[j+2]+2*(j+1)*c1[j+1];// Equation (5.9.2).
      cder2[j]=cder2[j+2]+2*(j+1)*c2[j+1];// Equation (5.9.2).
    }

  con=2.0/(XXX->get_xmax()-XXX->get_xmin());

  for (j=1;j<N_cheb;j++) //Normalize to the interval b-a.
    {
      cder1[j] *= con;
      cder2[j] *= con;
    }
  cder1[0] *= 1.0/(XXX->get_xmax()-XXX->get_xmin());
  cder2[0] *= 1.0/(XXX->get_xmax()-XXX->get_xmin());
}



void Derivateur::chder2(std::vector < std::vector<double> > const& c, std::vector < std::vector<double> > & cder, int iii)
//Given a,b,c[0..n-1], as output from routine chebft §5.8, and given n, the desired degree of approximation (length of c to be used), this routine returns the array cder[0..n-1], the Chebyshev coefficients of the derivative of the function whose coefficients are c.
{
  int j;
  double con;

  Chebychev *XXX = Chebychev::getInstance();
    
  cder[XXX->get_n()-1][iii]=0.0; //n-1 and n-2 are special cases.
  cder[XXX->get_n()-2][iii]=2*(XXX->get_n()-1)*c[XXX->get_n()-1][iii];
  for (j=XXX->get_n()-3;j>=0;j--)
    cder[j][iii]=cder[j+2][iii]+2*(j+1)*c[j+1][iii];// Equation (5.9.2).
  con=2.0/(XXX->get_xmax()-XXX->get_xmin());
  for (j=1;j<XXX->get_n();j++) //Normalize to the interval b-a.
    {
      cder[j][iii] *= con;
    }
  cder[0][iii] *= 1.0/(XXX->get_xmax()-XXX->get_xmin());
}




// First derivative calculation
double Derivateur::derivePrm(std::vector<double> const& fff, int iii, int nbPtRho)
{
  double fm2, fm1, f0, fp1, fp2;
  if(iii>=2&&iii<(nbPtRho-1))
    {
      fm2 = fff[iii-2];
      fm1 = fff[iii-1];
      f0 = fff[iii];
      fp1 = fff[iii+1];
      fp2 = fff[iii+2];
      return der1bulk(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==0)
    {
      fm2 = fff[0];
      fm1 = fff[1];
      f0 = fff[2];
      fp1 = fff[3];
      fp2 = fff[4];
      return der1FstPt(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==1)
    {
      fm2 = fff[0];
      fm1 = fff[1];
      f0 = fff[2];
      fp1 = fff[3];
      fp2 = fff[4];
      return der1ScdPt(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==(nbPtRho-1))
    {
      fm2 = fff[nbPtRho-4];
      fm1 = fff[nbPtRho-3];
      f0 = fff[nbPtRho-2];
      fp1 = fff[nbPtRho-1];
      fp2 = fff[nbPtRho];
      return der1LastPtm1(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==nbPtRho)
    {
      fm2 = fff[nbPtRho-4];
      fm1 = fff[nbPtRho-3];
      f0 = fff[nbPtRho-2];
      fp1 = fff[nbPtRho-1];
      fp2 = fff[nbPtRho];
      return der1LastPt(fm2,fm1,f0,fp1,fp2);
    }
  exit(0);
}




// Second derivative calculation
double Derivateur::deriveScd(std::vector<double> const& fff, int iii, int nbPtRho)
{
  double fm2, fm1, f0, fp1, fp2;
  if(iii>=2&&iii<(nbPtRho-1))
    {
      fm2 = fff[iii-2];
      fm1 = fff[iii-1];
      f0 = fff[iii];
      fp1 = fff[iii+1];
      fp2 = fff[iii+2];
      return der2bulk(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==0)
    {
      fm2 = fff[0];
      fm1 = fff[1];
      f0 = fff[2];
      fp1 = fff[3];
      fp2 = fff[4];
      return der2FstPt(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==1)
    {
      fm2 = fff[0];
      fm1 = fff[1];
      f0 = fff[2];
      fp1 = fff[3];
      fp2 = fff[4];
      return der2ScdPt(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==(nbPtRho-1))
    {
      fm2 = fff[nbPtRho-4];
      fm1 = fff[nbPtRho-3];
      f0 = fff[nbPtRho-2];
      fp1 = fff[nbPtRho-1];
      fp2 = fff[nbPtRho];
      return der2LastPtm1(fm2,fm1,f0,fp1,fp2);
    }
  else if(iii==nbPtRho)
    {
      fm2 = fff[nbPtRho-4];
      fm1 = fff[nbPtRho-3];
      f0 = fff[nbPtRho-2];
      fp1 = fff[nbPtRho-1];
      fp2 = fff[nbPtRho];
      return der2LastPt(fm2,fm1,f0,fp1,fp2);
    }
  exit(0);
}

