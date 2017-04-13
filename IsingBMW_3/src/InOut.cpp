#include "InOut.h"


InOut::InOut()
{

  // First three files opened
  fileInput.open("../input/data.ini", std::ios::in);
  log.open("../output/log");
}


void InOut::OpenFileRho0()
{
  fileRho0.open("../data/rho0", std::ios::app);
  fileRho0.precision(16);
}



void InOut::OpenSortiesGnuplot()
{
  fileW.open("../data/W", std::ios::app);
  fileW1.open("../data/W1", std::ios::app);
  fileW2.open("../data/W2", std::ios::app);
  fileV.open("../data/V", std::ios::app);

  fileDelta.open("../data/Delta", std::ios::app);
  fileDelta1.open("../data/Delta1", std::ios::app);
  fileDelta2.open("../data/Delta2", std::ios::app);

  fileI1.open("../data/I1", std::ios::app);
  fileI2.open("../data/I2", std::ios::app);
  fileI3.open("../data/I3", std::ios::app);
  fileI11.open("../data/I11", std::ios::app);
  fileJ3.open("../data/J3", std::ios::app);

  fileRegulateur.open("../data/Reg", std::ios::app);
  fileDerRegulateur.open("../data/DerReg", std::ios::app);
  filePropagateurQ.open("../data/Gq", std::ios::app);
  filePropagateurPQ.open("../data/Gppq", std::ios::app);

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

  fileRegulateur.precision(16);
  fileDerRegulateur.precision(16);
  filePropagateurQ.precision(16);
  filePropagateurPQ.precision(16);
  
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

  fileRegulateur.close();
  fileDerRegulateur.close();
  filePropagateurQ.close();
  filePropagateurPQ.close();
  
}




