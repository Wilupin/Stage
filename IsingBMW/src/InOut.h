#ifndef DEF_INOUT
#define DEF_INOUT

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>




class InOut
{

  

 public:
  
  // Fonctions qui en font un singleton
  static InOut  *getInstance()
  {
    if(NULL == _InOut)
      {
	_InOut = new InOut;
      }
    return _InOut;
  }

  static void kill()
  {
    if(NULL != _InOut)
      {
	delete _InOut;
	_InOut = NULL;
      }
  }

  
  // Fonctions de classe
  std::ifstream & get_fileInput() {return fileInput;};

  std::ofstream & get_log() {return log;};

  // Sorties gnuplot
  std::ofstream & get_fileW() {return fileW;};
  std::ofstream & get_fileW1() {return fileW1;};
  std::ofstream & get_fileW2() {return fileW2;};
  std::ofstream & get_fileV() {return fileV;};
  std::ofstream & get_fileDelta() {return fileDelta;};
  std::ofstream & get_fileDelta1() {return fileDelta1;};
  std::ofstream & get_fileDelta2(){return fileDelta2;};
  std::ofstream & get_fileI1() {return fileI1;};
  std::ofstream & get_fileI2() {return fileI2;};
  std::ofstream & get_fileI3() {return fileI3;};
  std::ofstream & get_fileI11() {return fileI11;};
  std::ofstream & get_fileJ3() {return fileJ3;};

  void OpenSortiesGnuplot(); 
  void CloseSortiesGnuplot();
  void OutEtaW0(int ttt, double eta, double w0, int irho0, double dt);

  void CloseFileInput() {fileInput.close();};
  void CloseLog() {log.close();};

  
 private:
  InOut();
  ~InOut(){};

  static InOut *_InOut;

  std::ifstream  fileInput;
  std::ofstream  log;
  
  // Sorties Gnuplot
  std::ofstream fileW;
  std::ofstream fileW1;
  std::ofstream fileW2;

  std::ofstream fileV;
  
  std::ofstream fileDelta;
  std::ofstream fileDelta1;
  std::ofstream fileDelta2;

  std::ofstream fileI1;
  std::ofstream fileI2;
  std::ofstream fileI3;
  std::ofstream fileI11;
  std::ofstream fileJ3;

  
};



#endif
