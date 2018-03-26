
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#ifndef MIXTUREPROFILEVI_H
#define MIXTUREPROFILEVI_H

#include <cmath>
#include "ProfileProcess.h"

// general superclass for all finite process mixtures on site-specific profiles
class MixtureProfileProcessVI: public virtual ProfileProcess	{

	public:

	MixtureProfileProcessVI() : profile(0) {}
	virtual ~MixtureProfileProcessVI(){}

	double* GetProfile(int site)	{
		return profile[alloc[site]];
	}

        double* GetStationaryDirweightVI(int site)    {

                return dirweightVI[alloc[site]];
        }

	int GetNcomponent() { return Ncomponent;}
	virtual int GetNDisplayedComponent()	{
		return Ncomponent;
	}
	int GetNOccupiedComponent() {
		int n = 0;
		for (int k=0; k<Ncomponent; k++)	{
			if (occupancy[k])	{
				n++;
			}
		}
		return n;
	}

	virtual int GetNmodeMax() {return GetNsite();}
	// virtual int GetNmodeMax() {return 100;}

	double GetMeanStationaryEntropy() {return GetStatEnt();}
	double GetSiteStationaryEntropy(int site) {return GetStatEnt(alloc[site]);}

	// summary statistic: mean entropy over all profiles
	double GetStatEnt();
	double GetStatEnt(int k);
	double GetCenterStatEnt();

	double GetMeanDirWeight();

	void RenormalizeProfiles();

	// generic Move function
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1) = 0;
	double MoveDirWeights(double tuning, int nrep);

	protected:

	virtual void UpdateModeProfileSuffStat() = 0;

	// implements a pure virtual defined in ProfileProcess
	double ProfileSuffStatLogProb();

	// suffstat lnL of site <site> when allocated to component <cat>
	virtual double LogStatProb(int site, int cat) = 0;

	// the component suff stat log prob is yet to be implemented in subclasses
	virtual double ProfileSuffStatLogProb(int cat) = 0;

	// called at the beginning and end of the run (see PhyloProcess)
	virtual void Create(int innsite, int indim);
	virtual void Delete();

	// in certain models,
	// the matrix associated to each component should be created on the fly
	// all other component-specific variables are static (see above, NmodeMax)
	virtual void CreateComponent(int k) = 0;
	// virtual void CreateComponent(int k, double* instat) = 0;
	virtual void DeleteComponent(int k) = 0;
	virtual void UpdateComponent(int k) = 0;

	void UpdateComponents()	{
		for (int k=0; k<GetNcomponent(); k++)	{
			UpdateComponent(k);
		}
	}

	virtual double GetAllocEntropy()	{
		double total = 0;
		UpdateOccupancyNumbers();
		for (int k=0; k<GetNcomponent(); k++)	{
			double tmp = ((double) occupancy[k]) / GetNsite();
			if (tmp)	{
				total -= tmp * log(tmp);
			}
		}
		return total;
	}

	virtual void AddSite(int site, int cat)	{
		alloc[site] = cat;
		occupancy[cat]++;
	}
	virtual void RemoveSite(int site, int cat)	{
		occupancy[cat]--;
	}

	// sample all aspects of the mixture (number of components, composition) from the prior
	void SampleProfile();

	virtual void SampleHyper() = 0;
	virtual void SampleAlloc() = 0;
	virtual void SampleStat();
	void SampleStat(int cat);
	void SampleStat(double* stat, double statmin = 0);

	void UpdateOccupancyNumbers();
	double ResampleEmptyProfiles();

	double LogProfilePrior();

	virtual double LogHyperPrior() = 0;
	virtual double LogAllocPrior()	{
		return 0;
	}

	// dirichlet prior ~ Dirichlet (dirweight[0]... dirweight[GetDim()-1])
	virtual double LogStatPrior();

	virtual double LogStatPrior(int cat);

	virtual void SwapComponents(int cat1, int cat2);
        // virtual double GetPCAT(int cat) = 0;      
 
        // double* PCAT; 
        double* LearningDirweight; 
        double** dirweightVI;
        double** dirweighthat;
        double** dirweight_grad;
	double** profile;
	double* allocprofile;
	double* dirweight;
	int* alloc;
	int* occupancy;
	int Ncomponent;
	double* logstatprior;
	double* profilesuffstatlogprob;
};

#endif

