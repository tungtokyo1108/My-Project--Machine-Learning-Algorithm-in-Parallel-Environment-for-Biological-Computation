
/********************


**********************/


#ifndef POISSONSUB_H
#define POISSONSUB_H

#include "SubstitutionProcess.h"
#include "PoissonProfileProcess.h"
#include "RateProcessVI.h"
#include <gsl/gsl_sf_psi.h>
#include "Random.h"

class PoissonSubstitutionProcess : public virtual SubstitutionProcess, public virtual PoissonProfileProcess, public virtual RateProcessVI {

	public:

	PoissonSubstitutionProcess() : zipstat(0) {}
	virtual ~PoissonSubstitutionProcess() {}


	const double* GetStationary(int site) {return zipstat[site];}
	int GetNstate(int site) {return GetZipSize(site);}
        // virtual void GetRandomStateFromZip(int* stateup, int* statedown, double time);

	protected:

	// CPU Level 3: implementations of likelihood propagation and substitution mapping methods
	void Propagate(double** from, double** to, double time, bool condalloc = false);
	BranchSitePath** SamplePaths(int* stateup, int* statedown, double time);
        BranchSitePath** SamplePathsVI(int* stateup, int* statedown, double timealpha, double timebeta);
	BranchSitePath** SampleRootPaths(int* rootstate);
         

	void SimuPropagate(int* stateup, int* statedown, double time);
	void SimuPropagateZip(int* stateup, int* statedown, double time);

	// CPU Level 1: gathering sufficient statistics from substitution mappings
	void AddSiteRateSuffStat(int* siteratesuffstatcount, BranchSitePath** patharray);
	void AddBranchLengthSuffStat(int& count, BranchSitePath** patharray);
	void AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root);
        // void AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, double numbertypesubstitution, bool root);
        // int typeaminoacid;

	virtual void Delete() {
		DeleteZip();
		SubstitutionProcess::Delete();
	}

	void CreateZip();
	void DeleteZip();
	void UpdateZip();
	void UpdateZip(int site);
	void UnzipBranchSitePath(BranchSitePath** patharray, int* nodestateup, int* nodestatedown);
	int GetRandomStateFromZip(int site, int state);

	// will be implemented in PoissonPhyloProcess
	virtual int GetZipSize(int site) = 0;
	virtual int GetOrbitSize(int site) = 0;
	virtual int GetStateFromZip(int site, int state) = 0;
	virtual bool InOrbit(int site, int state) = 0;

	void ChooseTrueStates(BranchSitePath** patharray, int* nodestateup, int* nodestatedown, bool root);
	void ChooseRootTrueStates(int* nodestate);

	private:

	double** zipstat;
};

#endif

