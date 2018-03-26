
/********************


**********************/

#ifndef POISSONMIXTUREPROFILEVI_H
#define POISSONMIXTUREPROFILEVI_H

#include "PoissonProfileProcess.h"
#include "MixtureProfileProcessVI.h"
#include <gsl/gsl_sf_psi.h>

// superclass for Poisson (F81) implementations
class PoissonMixtureProfileProcessVI: public virtual PoissonProfileProcess, public virtual MixtureProfileProcessVI	{

	public:

	PoissonMixtureProfileProcessVI() : profilesuffstatcount(0) {}
	virtual ~PoissonMixtureProfileProcessVI() {}

        // virtual double GetPCAT(int cat) = 0;
        // virtual double GetTotalSitePCAT() = 0;
        double Getprofilesuffstatcount(int cat);
        // double Getdirweighthat(int site) {return dirweighthat[site];}
        // double** dirweightVI;
        // double** dirweighthat;
        virtual double GetProfileParameter(int site, int cat);
        double TotDirWeightVI(int cat, int dim);
        double PCATDirWeight();
        double PCATDirWeight(int cat);
        double ELBODirWeight();
        double ELBODirWeight(int cat);
        double ELBO_DirWeight;        

	protected:

	virtual void Create(int innsite, int indim);
	virtual void Delete();

	virtual void CreateComponent(int k) {
		occupancy[k] = 0;
		/*
		int* catnsub = profilesuffstatcount[k];
		for (int i=0; i<GetDim(); i++)	{
			catnsub[i] = 0;
		}
		*/
		SampleStat(k);
	}
	virtual void DeleteComponent(int k) {
	}
	virtual void UpdateComponent(int k) {}

	// posterior
	// collects sufficient statistics across sites, pools them componentwise
	void UpdateModeProfileSuffStat();

	// virtual void CreateComponent(int k)	{}

	// suffstat lnL of all sites allocated to component cat
	double ProfileSuffStatLogProb(int cat);

	// difference between:
	// suffstat lnL of all sites allocated to component cat when site <site> is among them, and
	// suffstat lnL of all sites allocated to component cat when site <site> is not among them
	double DiffLogSampling(int cat, int site);
	virtual double LogStatProb(int site, int cat);
        double EstimateDirWeight();
	double MoveDirWeightVI();

	double MoveProfile();
	double MoveProfile(int cat);

	void SwapComponents(int cat1, int cat2);
	void AddSite(int site, int cat);
	void RemoveSite(int site, int cat);

	double GetNormRate(int k)	{

		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			for (int j=i+1; j<GetDim(); j++)	{
				tot += profile[k][i] * profile[k][j];
			}
		}
		return 2*tot;
	}

	virtual double GetNormalizationFactor()	{
		UpdateOccupancyNumbers();
		double norm = 0;
		int tot = 0;
		for (int k=0; k<GetNcomponent(); k++)	{
			if (occupancy[k])	{
				double tmp = GetNormRate(k);
				norm += (occupancy[k] + 1) * tmp;
				tot += occupancy[k] + 1;
			}
		}
		/*
		if (tot != GetNsite() + GetNcomponent())	{
			cerr << "error in norm factor\n";
			cerr << tot << '\t' << GetNsite() << '\t' << GetNcomponent() << '\n';
			exit(1);
		}
		*/
		norm /= tot;
		return norm;
	}
        
	// private:
        /*virtual double GetPCAT(int cat) = 0;
        double Getprofilesuffstatcount(int cat);
        // double Getdirweighthat(int site) {return dirweighthat[site];}
        double** dirweightVI;
        double** dirweighthat;
        virtual double GetProfileParameter(int site, int cat);*/

	int** profilesuffstatcount;
};

#endif

