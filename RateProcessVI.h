
/********************

**********************/


#ifndef RATEVI_H
#define RATEVI_H

#include <iostream>
using namespace std;

#include "Chrono.h"

class RateProcessVI {

	public:

	RateProcessVI() : nsite(0) {}
	virtual ~RateProcessVI() {}

	virtual string GetVersion() = 0;

	int GetNsite() {return nsite;}

	virtual int GetNrate(int site) = 0;
	virtual double GetRate(int site) = 0;
	virtual double GetRateWeight(int site, int cat) = 0;
	double GetMeanRateVI();
        double GetMeanRate();
        
	virtual double GetPriorMeanRate() = 0;
	virtual double GetAlpha() {return 1;}

	virtual void ActivateSumOverRateAllocations() = 0;
	virtual void InactivateSumOverRateAllocations(int* ratealloc) = 0;
	bool SumOverRateAllocations() {return sumflag;}

	virtual double LogRatePrior() = 0;
	virtual void SampleRate() = 0;

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;

	// protected:

	// abstract classes will be implemented in phyloprocess
	virtual void GlobalUpdateSiteRateSuffStat() = 0;
	virtual void SlaveUpdateSiteRateSuffStat() = 0;

	virtual void UpdateSiteRateSuffStat() = 0;
	virtual double GetSiteRateSuffStatBeta(int site) = 0;
	virtual int GetSiteRateSuffStatCount(int site) = 0;
        // virtual int GetSiteRateSuffStatCountELBO(int site) = 0;

        /*virtual void TotalParameterLength() = 0;
        virtual double GetTotalParaLength() = 0;*/

	void Create(int innsite)	{
		nsite = innsite;
	}
	void Delete() {}

	virtual int GetNprocs() = 0;
	virtual int GetMyid() = 0;
	virtual int GetSiteMin() = 0;
	virtual int GetSiteMax() = 0;

	bool sumflag;
	int nsite;
        double MeanRateVI;

        // double* RateAlphaVI;
        // double* RateBetaVI;
        /*virtual void TotalParameterRate() = 0;
        double TotalParaRate;
        virtual double GetTotalParaRate() = 0;*/
        virtual double GetRateAlphaVI(int site) = 0;
        virtual double GetRateBetaVI(int site) = 0;
        virtual double GlobalGetMeanRateVI() = 0;
        virtual void SlaveGetMeanRateVI() = 0;
         
        virtual void GlobalGetRateAlphaVI() = 0;
        virtual void SlaveGetRateAlphaVI() = 0;
        virtual void GlobalGetRateBetaVI() = 0;
        virtual void SlaveGetRateBetaVI() = 0;

	Chrono chronorate;
};

#endif

