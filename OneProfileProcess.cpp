
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#include "OneProfileProcess.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* OneProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void OneProfileProcess::Create(int innsite, int indim)	{
	if (! profile)	{
		ProfileProcess::Create(innsite,indim);
		profile = new double[GetDim()];
	}
}

void OneProfileProcess::Delete()	{
	if (profile)	{
		delete[] profile;
		profile = 0;
		ProfileProcess::Delete();
	}
}

double OneProfileProcess::GetStatEnt()	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		total -= profile[i] * log(profile[i]);
	}
	return  total;
}

void OneProfileProcess::SampleProfile()	{
	SampleStat();
}

void OneProfileProcess::SampleStat()	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[k] = rnd::GetRandom().sGamma(1.0);
		total += profile[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[k] /= total;
	}
	total = 0;
	for (int k=0; k<GetDim(); k++)	{
		if (profile[k] < stateps)	{
			profile[k] = stateps;
		}
		total += profile[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[k] /= total;
	}
	// UpdateComponent(i);
}

double OneProfileProcess::LogProfilePrior()	{
	return 0;
}


