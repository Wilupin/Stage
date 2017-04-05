#include "InOut.h"


InOut::InOut()
{

  // Two first files opened
  fileInput.open("../input/data.ini", std::ios::in);
  log.open("../output/log");
}


void InOut::OpenSortiesGnuplot()
{
  fileW.open("../data/W", std::ios::app);
  fileW1.open("../data/W", std::ios::app);
  fileW2.open("../data/W", std::ios::app);
  fileV.open("../data/W", std::ios::app);

  fileDelta.open("../data/W", std::ios::app);
  fileDelta1.open("../data/W", std::ios::app);
  fileDelta2.open("../data/W", std::ios::app);

  fileI1.open("../data/W", std::ios::app);
  fileI2.open("../data/W", std::ios::app);
  fileI3.open("../data/W", std::ios::app);
  fileI11.open("../data/W", std::ios::app);
  fileJ3.open("../data/W", std::ios::app);


  fileW.precision(16);
  fileW1.precision(16);
  fileW2.precision(16);
  fileV.precision(16);

  fileDelta.precision(16);
  fileDelta1.precision(16);
  fileDelta2.precision(16);

  fileI1.precision(16);
  fileI2.precision(16);
  fileI3.precision(16);
  fileI11.precision(16);
  fileJ3.precision(16);
  
}



void InOut::CloseSortiesGnuplot()
{
  fileW.close();
  fileW1.close();
  fileW2.close();
  fileV.close();

  fileDelta.close();
  fileDelta1.close();
  fileDelta2.close();

  fileI1.close();
  fileI2.close();
  fileI3.close();
  fileI11.close();
  fileJ3.close();
  
}



void InOut::OutEtaW0(int ttt, double eta, double w0, int irho0, double dt)
{
    std::ofstream w0etaout("../data/w0eta",std::ios::app);
    w0etaout.precision(16);
    w0etaout<<ttt*dt<<"\t"<<eta<<"\t"<<w0<<"\t"<<irho0<<std::endl;
    w0etaout.close();
}
