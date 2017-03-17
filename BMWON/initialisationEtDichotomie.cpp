#include "initialisationEtDichotomie.h"
#include "function.h"

void ReadAndInitialise(GammaAndDerivatives & gamma, GaussLegendre & gausslegvar, DichoParameter & dichoprm)
{
    system("rm -f data/*");
    ofstream log("log");

    ifstream dataIni("data.ini");
    if(dataIni==NULL) log << "erreur : fichier d'entree introuvable !!\n" << endl, log.close(), exit(0);
    int NUM_THREADS = 0;

    string chaine="";
    getline(dataIni,chaine);
    dataIni>>gausslegvar.d>>gamma.n>>gausslegvar.alpha>>dichoprm.u0>>dichoprm.rmin>>dichoprm.rmax>>dichoprm.rht>>dichoprm.rlt>>dichoprm.OnOff>>gamma.tmax>>gamma.ompOn>>NUM_THREADS;
    dataIni.close();

    log<<"Programme BMW O(n = "<<gamma.n<<") en d = "<<gausslegvar.d<<" (v03.12.15)"<<endl;
    log<<"nppp = "<<nppp<<";\t nrho = "<<nrho<<";\t ngauleg = "<<ngausslegendre<<endl;
    log<<"pmax = "<<p_max<<";\t rhomax = "<<rho_max<<";\t qInt = "<< qInt<<endl;
    log<<"alpha = "<<gausslegvar.alpha<<";\t u0 = "<<dichoprm.u0<<endl;
    gamma.d = gausslegvar.d;
    gamma.pPrescription = nppp-1, gamma.rhoPrescription = 5;
    ChebRoot(gamma.Pcheb,"data/ppp",npcheb,npchebev,p_min,p_max);
    gaulegWeightAbscissas(gausslegvar);

    log<<"Prescription Ya(" <<gamma.Pcheb.xcheb[gamma.pPrescription]<<" , "<< gamma.rhoPrescription*drho <<" ) = 0.0"<<endl;

    initialiseGamma(gamma);

    dichoprm.tht = 0, dichoprm.tlt = 0, dichoprm.nrun = 0, dichoprm.avc = 0;
    dichoprm.tttavc = 80000;
    dichoprm.seuil = 1e-8;

    if(gamma.ompOn)
    {
        omp_set_num_threads(NUM_THREADS);
        log<<"On parallélise sur "<<NUM_THREADS<<" processeurs.\n"<<endl;
    }
    log.close();
}

void ChebRoot(chebvar& XXX, char OutFile[], int ncheb, int nevcheb, double x_min, double x_max)
{
    XXX.n = ncheb; XXX.nev = nevcheb; XXX.xmax=x_max; XXX.xmin = x_min;
    ofstream OUTRoot(OutFile);
    OUTRoot.precision(16);
    int kkk;
    double x, bma = 0.5*(XXX.xmax-XXX.xmin), bpa = 0.5*(XXX.xmax+XXX.xmin);
    for(kkk=0;kkk<XXX.n;kkk++)
    {
        x = cos(M_PI*(kkk+0.5)/XXX.n);
        XXX.xcheb.push_back(x*bma+bpa);
        OUTRoot<<XXX.xcheb[kkk]<<endl;
    }
    OUTRoot.close();
}

void gaulegWeightAbscissas(GaussLegendre& gausslegvar)
{
    int m,i,j;
    double z1,z,pp,p3,p2,p1; //High precision is a good idea for this routine.

    gausslegvar.x.resize(gausslegvar.ngl+2);
    gausslegvar.w.resize(gausslegvar.ngl+2);

    m=(gausslegvar.ngl+1)/2;                     //The roots are symmetric in the interval, so we only have to find half of them
    for(i=1;i<=m;i++)
    //Loop over the desired roots
    {
        z = cos(M_PI*(i-0.25)/(gausslegvar.ngl+0.5)); //Starting with the above approximation the ith root, we enter the main loop of refinement by Newtow's method.
        do
        {
            p1=1.0;
            p2=0.0;
            for(j=1;j<=gausslegvar.ngl;j++)//Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
            {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            //p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2, the polynomial of one lower order.
            pp=gausslegvar.ngl*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;     //Newthon's method.
        }while(fabs(z-z1)>EPS);
        gausslegvar.x[i]=-z;                       //Scale the root to the desired interval, and put in its symmetric counterpart.
        gausslegvar.x[gausslegvar.ngl+1-i]=z;
        gausslegvar.w[i]=2./((1.0-z*z)*pp*pp);      //Compute the weight and its symmetric counterpart
        gausslegvar.w[gausslegvar.ngl+1-i]=gausslegvar.w[i];
    }
    gausslegvar.norme = normalisationIntegrale(gausslegvar.d);
}

double normalisationIntegrale(double d)
{
    return qInt*0.5*qInt*d*tgamma(d/2.)/tgamma((d-1)/2.)/sqrt(M_PI);
}

void initialiseGamma(GammaAndDerivatives& gamma)
{
    gamma.YYYa_p_rho.resize(nppp),gamma.YYYa1_p_rho.resize(nppp),gamma.YYYa2_p_rho.resize(nppp),gamma.YYYb_p_rho.resize(nppp),gamma.YYYb1_p_rho.resize(nppp),gamma.YYYb2_p_rho.resize(nppp);
    gamma.JJJ3LL_p_rho.resize(nppp),gamma.JJJ3TT_p_rho.resize(nppp),gamma.JJJ3LT_p_rho.resize(nppp),gamma.JJJ3TL_p_rho.resize(nppp);
    gamma.chebcoef_YYYa1_rho_p.resize(nrho+1),gamma.chebcoef_YYYa_rho_p.resize(nrho+1), gamma.chebcoef_YYYap1_rho_p.resize(nrho+1),gamma.chebcoef_YYYb1_rho_p.resize(nrho+1),gamma.chebcoef_YYYb_rho_p.resize(nrho+1), gamma.chebcoef_YYYbp1_rho_p.resize(nrho+1);

    gamma.III2LL_rho.resize(nrho+1,0.0),gamma.III2TT_rho.resize(nrho+1,0.0),gamma.III3LL_rho.resize(nrho+1,0.0),gamma.III3TT_rho.resize(nrho+1,0.0),gamma.III3LT_rho.resize(nrho+1,0.0),gamma.III3TL_rho.resize(nrho+1,0.0),gamma.IIIA_rho.resize(nrho+1,0.0),gamma.WWW_rho.resize(nrho+1,0.0),gamma.WWW1_rho.resize(nrho+1,0),gamma.WWW2_rho.resize(nrho+1,0),gamma.WWWI1_rho.resize(nrho+1,0);

    for(int ipp = 0; ipp < nppp; ipp++) gamma.YYYa_p_rho[ipp].resize(nrho+1,0.0),gamma.YYYa1_p_rho[ipp].resize(nrho+1,0.0),gamma.YYYa2_p_rho[ipp].resize(nrho+1,0.0),gamma.YYYb_p_rho[ipp].resize(nrho+1,0.0),gamma.YYYb1_p_rho[ipp].resize(nrho+1,0.0),gamma.YYYb2_p_rho[ipp].resize(nrho+1,0.0),gamma.JJJ3LL_p_rho[ipp].resize(nrho+1,0.0),gamma.JJJ3TT_p_rho[ipp].resize(nrho+1,0.0),gamma.JJJ3LT_p_rho[ipp].resize(nrho+1,0.0),gamma.JJJ3TL_p_rho[ipp].resize(nrho+1,0.0);

    for(int iii = 0; iii <= nrho; iii++)  gamma.chebcoef_YYYa1_rho_p[iii].resize(npcheb,0.0),gamma.chebcoef_YYYa_rho_p[iii].resize(npcheb,0.0), gamma.chebcoef_YYYap1_rho_p[iii].resize(npcheb,0.0),gamma.chebcoef_YYYb1_rho_p[iii].resize(npcheb,0.0),gamma.chebcoef_YYYb_rho_p[iii].resize(npcheb,0.0), gamma.chebcoef_YYYbp1_rho_p[iii].resize(npcheb,0.0);
}

int initialising(GammaAndDerivatives& gamma, DichoParameter & dichoprm)
{
    dichoprm.nrun++;

    if(dichoprm.avc)
    {
        dichoprm.r0 = 0.5*(dichoprm.rmax + dichoprm.rmin);
        gamma.eta = dichoprm.eta_avc;
        gamma.Zk = dichoprm.Zk_avc;
        gamma.kk = exp(dichoprm.tttavc*dtt);
        for(int iii = 0; iii <= nrho; iii++)
        {
            gamma.WWW_rho[iii] = dichoprm.WWW_rho_avc[iii]+dichoprm.r0-dichoprm.w0_avc;
            for(int ipp = 0; ipp < nppp; ipp++) gamma.YYYa_p_rho[ipp][iii]=dichoprm.YYYa_p_rho_avc[ipp][iii],gamma.YYYb_p_rho[ipp][iii]=dichoprm.YYYb_p_rho_avc[ipp][iii];
        }
//        gamma.irho_0 = set_irho_0(gamma.WWW_rho);
//        gamma.rho_0 = ZeroW(gamma);
        return dichoprm.tttavc;
    }
    else
    {
        dichoprm.r0 = 0.5*(dichoprm.rmax + dichoprm.rmin);
        gamma.eta = 0.0;
        gamma.Zk=1.0;
        gamma.kk=1.0;
        for(int iii = 0; iii <= nrho; iii++)
        {
            gamma.WWW_rho[iii] = dichoprm.r0 + dichoprm.u0*iii*drho;
            for(int ipp = 0; ipp < nppp; ipp++) gamma.YYYa_p_rho[ipp][iii]=0.0,gamma.YYYb_p_rho[ipp][iii]=0.0;
        }
//        gamma.irho_0=set_irho_0(gamma.WWW_rho);
//        gamma.rho_0 = ZeroW(gamma);
        return 0;
    }
}

int dichotomie(double w0, DichoParameter & dichoprm, int ttt, ofstream & log)
{
    ofstream w0etaout("data/w0eta",ios::app);
    w0etaout<<endl;
    w0etaout.close();

    if(w0<dichoprm.rlt)
    {
        if(ttt>dichoprm.tlt)
        {
            log<<"Dichotomie à t = "<<ttt*dtt<<endl<<"Phase basse temperature rmin -> r0 "<<endl;
            dichoprm.tlt = ttt;
            dichoprm.rmin = dichoprm.r0;
        }
        else
        {
            log<<"La dichotomie n'avance plus, fin du programme!"<<endl;
            log<<"Paramètre dichotomie fin :\n rmin = " << dichoprm.rmin<<",\n rmax = " << dichoprm.rmax<<",\n u0 = " << dichoprm.u0<<",\n tht = " << dichoprm.tht<<",\n tlt = " << dichoprm.tlt<<endl;
            log.close();
            exit(0);
        }
    }
    else
    {
        if(ttt>dichoprm.tht)
        {
            log<<"Dichotomie à t = "<<ttt*dtt<<endl<<"Phase haute temperature rmax -> r0"<<endl;
            dichoprm.tht = ttt;
            dichoprm.rmax = dichoprm.r0;
        }
        else
        {
            log<<"La dichotomie n'avance plus, fin du programme!"<<endl;
            log<<"Paramètre dichotomie fin :\n rmin = " << dichoprm.rmin<<",\n rmax = " << dichoprm.rmax<<",\n u0 = " << dichoprm.u0<<",\n tht = " << dichoprm.tht<<",\n tlt = " << dichoprm.tlt<<endl;
            log.close();
            exit(0);
        }
    }
    return -1;
}


void sortieDiversesGnuPlot(GammaAndDerivatives const& gamma, int ttt, int nrun)
{
    ofstream outfffa("data/fffa",ios::app),outfffa1("data/fffa1",ios::app),outfffa2("data/fffa2",ios::app),outfffap1("data/fffap1",ios::app),outfffb("data/fffb",ios::app),outfffb1("data/fffb1",ios::app),outfffb2("data/fffb2",ios::app),outfffbp1("data/fffbp1",ios::app),outWWW("data/W",ios::app),outWWW1("data/W1",ios::app),outWWW2("data/W2",ios::app), outdI1("data/dI1",ios::app), outI2ll("data/I2ll",ios::app), outI2tt("data/I2tt",ios::app), outI3ll("data/I3ll",ios::app), outI3tt("data/I3tt",ios::app), outI3lt("data/I3lt",ios::app), outI3tl("data/I3tl",ios::app), outJ3ll("data/J3ll",ios::app), outJ3tt("data/J3tt",ios::app), outJ3lt("data/J3lt",ios::app), outJ3tl("data/J3tl",ios::app), outIa("data/Ia");
    outfffa.precision(16),outfffa1.precision(16),outfffa2.precision(16),outfffap1.precision(16),outfffb.precision(16),outfffb1.precision(16),outfffb2.precision(16),outfffbp1.precision(16),outWWW.precision(16),outWWW1.precision(16),outWWW2.precision(16),outdI1.precision(16),outI2ll.precision(16),outI2tt.precision(16),outI3ll.precision(16),outI3tt.precision(16),outI3lt.precision(16),outI3tl.precision(16),outJ3ll.precision(16),outJ3tt.precision(16),outJ3lt.precision(16),outJ3tl.precision(16),outIa.precision(16);

    for(int iii=0; iii <= nrho; iii++)
    {
        outdI1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.WWWI1_rho[iii]<<endl;
        outIa<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.IIIA_rho[iii]<<endl;
        outI2ll<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.III2LL_rho[iii]<<endl;
        outI2tt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.III2TT_rho[iii]<<endl;
        outI3ll<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.III3LL_rho[iii]<<endl;
        outI3tt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.III3TT_rho[iii]<<endl;
        outI3lt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.III3LT_rho[iii]<<endl;
        outI3tl<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.III3TL_rho[iii]<<endl;
        outWWW <<nrun<<"\t"<<ttt<<"\t"<< iii*drho<<"\t"<<gamma.WWW_rho[iii]<<endl;
        outWWW1 <<nrun<<"\t"<<ttt<<"\t"<< iii*drho<<"\t"<<gamma.WWW1_rho[iii]<<endl;
        outWWW2 <<nrun<<"\t"<<ttt<<"\t"<< iii*drho<<"\t"<<gamma.WWW2_rho[iii]<<endl;

        for(int ipp=0 ; ipp < nppp ; ipp++)
        {
            outfffa<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.YYYa_p_rho[ipp][iii]<<endl;
            outfffa1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.YYYa1_p_rho[ipp][iii]<<endl;
            outfffa2<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.YYYa2_p_rho[ipp][iii]<<endl;
            outfffap1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<chebev(gamma.chebcoef_YYYap1_rho_p[iii], gamma.Pcheb,gamma.Pcheb.xcheb[ipp])<<endl;
            outfffb<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.YYYb_p_rho[ipp][iii]<<endl;
            outfffb1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.YYYb1_p_rho[ipp][iii]<<endl;
            outfffb2<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.YYYb2_p_rho[ipp][iii]<<endl;
            outfffbp1<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<chebev(gamma.chebcoef_YYYbp1_rho_p[iii], gamma.Pcheb,gamma.Pcheb.xcheb[ipp])<<endl;
            outJ3ll<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.JJJ3LL_p_rho[ipp][iii]<<endl;
            outJ3tt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.JJJ3TT_p_rho[ipp][iii]<<endl;
            outJ3lt<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.JJJ3LT_p_rho[ipp][iii]<<endl;
            outJ3tl<<nrun<<"\t"<<ttt<<"\t"<<iii*drho<<"\t"<<gamma.Pcheb.xcheb[ipp]<<"\t"<<gamma.JJJ3TL_p_rho[ipp][iii]<<endl;
        }
        outfffa<<endl;
        outfffa1<<endl;
        outfffa2<<endl;
        outfffap1<<endl;
        outfffb<<endl;
        outfffb1<<endl;
        outfffb2<<endl;
        outfffbp1<<endl;
        outJ3ll<<endl;
        outJ3tt<<endl;
        outJ3lt<<endl;
        outJ3tl<<endl;
    }

    outfffa.close(),outfffa1.close(),outfffa2.close(),outfffap1.close(),outfffb.close(),outfffb1.close(),outfffb2.close(),outfffbp1.close(),outWWW.close(),outWWW1.close(),outWWW2.close(),outdI1.close(),outI2ll.close(),outI2tt.close(),outI3ll.close(),outI3tt.close(),outI3lt.close(),outI3tl.close(),outIa.close(),outJ3ll.close(),outJ3tt.close(),outJ3lt.close(),outJ3tl.close();
}

