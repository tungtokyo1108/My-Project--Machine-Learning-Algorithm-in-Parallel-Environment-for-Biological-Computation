
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#ifndef ONEPROFILE_H
#define ONEPROFILE_H

#include <cmath>
#include "ProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class OneProfileProcess: public virtual ProfileProcess	{

	public:

	OneProfileProcess() : profile(0) {}
	virtual ~OneProfileProcess(){}

	double* GetProfile(int site)	{
		return profile;
	}

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt();}

	int GetNOccupiedComponent()  {return 1;}

	double GetStatEnt();

	// generic Move function
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1) = 0;

	protected:

	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create(int innsite, int indim);
	virtual void Delete();

	void UpdateProfile()	{
	}


	// sample all aspects of the mixture (number of components, composition) from the prior
	void SampleProfile();
	void SampleStat();

	double LogProfilePrior();

	double* profile;
};

#endif

