#include "Parametres.h"


Parametres::Parametres()
{

  std::ifstream fichierInput("../input/data.ini", std::ios::in);
  
  
  if(fichierInput)
    {
      std::string ligne;
      std::string param;

      while(getline(fichierInput, ligne))
	{
	  
	  size_t pos = ligne.find("=");
	  
	  if(pos != -1)
	    {
	      
	      param = trim(ligne.substr(0,pos));

	      // Parametres generaux
	      if (param == "alpha")
		alpha = std::stod(trim(ligne.substr(pos+1)));
	      if (param == "mu")
		mu = std::stod(trim(ligne.substr(pos+1)));

	      // Parametre temperature
	      if (param == "beta")
		beta = std::stod(trim(ligne.substr(pos+1)));

	      // parmetres sur le temps
	      if (param == "tMax")
		tMax = std::stod(trim(ligne.substr(pos+1)));
	      if (param == "tMin")
		tMin = std::stod(trim(ligne.substr(pos+1)));
	      if (param == "dt")
		dt = std::stod(trim(ligne.substr(pos+1)));

	      // parametres de paralellisation
	      if (param == "ompOn")
		ompOn = std::stoi(trim(ligne.substr(pos+1)));
	      if (param == "numberOfThreads")
		numberOfThreads = std::stoi(trim(ligne.substr(pos+1)));

	      // Grille sur rho
	      if (param == "rhoMin")
		rhoMin = std::stod(trim(ligne.substr(pos+1)));
	      if (param == "rhoMax")
		rhoMax = std::stod(trim(ligne.substr(pos+1)));
	      if (param == "drho")
		drho = std::stod(trim(ligne.substr(pos+1)));

	      // Parametres pour Gauss Legendre
	      if (param == "qMax")
		qMax = std::stod(trim(ligne.substr(pos+1)));
	      if (param ==  "npGL")
		npGL = std::stoi(trim(ligne.substr(pos+1)));

	      // Parametres pour Chebychev
	      if (param == "npCheb")
		npCheb = std::stoi(trim(ligne.substr(pos+1)));
	      if (param == "pMax")
		pMax = std::stod(trim(ligne.substr(pos+1)));
	      if (param == "pMin")
		pMin = std::stod(trim(ligne.substr(pos+1)));

	      // Parametres d'affichage
	      if (param == "itSaveInLog")
		itSaveInLog = std::stoi(trim(ligne.substr(pos+1)));
	      if (param == "itSaveInData")
		itSaveInData = std::stoi(trim(ligne.substr(pos+1)));
	      if (param == "itSaveEvthgInData")
		itSaveEvthgInData = std::stoi(trim(ligne.substr(pos+1)));
	    }
	 	  
	}
    }

  

  // Calculs des coeffcients restants
  tMaxInt   = int(-(tMax-tMin)/dt);  
  npChebev  = int(3*npCheb/4);
  nrho      = int((rhoMax-rhoMin)/drho);

}





std::string Parametres::trim(std::string const& str)
{
  size_t first = str.find_first_not_of(' ');
  size_t last  = str.find_last_not_of(' ');

  return str.substr(first, (last-first+1));
}
