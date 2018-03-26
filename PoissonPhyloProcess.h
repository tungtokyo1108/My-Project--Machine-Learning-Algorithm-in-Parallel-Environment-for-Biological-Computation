
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#ifndef POISSONPHYLO_H
#define POISSONPHYLO_H

#include "PhyloProcess.h"
#include "PoissonSubstitutionProcess.h"

class PoissonPhyloProcess : public virtual PhyloProcess, public virtual PoissonSubstitutionProcess	{

	public:

	PoissonPhyloProcess() : siteprofilesuffstatcount(0), allocsiteprofilesuffstatcount(0), zipdata(0) {}
	virtual ~PoissonPhyloProcess() {}

	void Unfold()	{
		UpdateZip();
		PhyloProcess::Unfold();
	}

	// protected:

	// true data here !
	virtual void Create(Tree* intree, SequenceAlignment* indata);
	virtual void Delete();

	// in fact, same object as GetData, but now with its true type
	ZippedSequenceAlignment* GetZipData()	{
		if (! zipdata)	{
			cerr << "null zip\n";
			exit(1);
		}
		return zipdata;
	}

	virtual int GetNstate() {return truedata->GetNstate();}
	int GetZipSize(int site) {return GetZipData()->GetZipSize(site);}
	int GetOrbitSize(int site) {return GetZipData()->GetOrbitSize(site);}
	int GetStateFromZip(int site, int state) {return GetZipData()->GetStateFromZip(site,state);}
	bool InOrbit(int site, int state) {return GetZipData()->InOrbit(site,state);}
	
	void CreateSuffStat();
	void DeleteSuffStat();

	void UpdateSiteRateSuffStat();
	void UpdateSiteProfileSuffStat();
	void UpdateBranchLengthSuffStat();

	void GlobalUpdateSiteProfileSuffStat();
	void SlaveUpdateSiteProfileSuffStat();

	const int* GetSiteProfileSuffStatCount(int site) {return siteprofilesuffstatcount[site];}

	void SetDataFromLeaves()	{
		SampleTrueNodeStates(GetRoot());
		PhyloProcess::SetDataFromLeaves();
	}

	// virtual void RecursiveSimulateForward(const Link* from);

	void SampleTrueNodeStates(const Link* from);

	virtual double GetObservedCompositionalHeterogeneity(double* taxstat, double& meandist)	{
		return truedata->CompositionalHeterogeneity(taxstat,0,meandist);
	}

	void RecursiveUnzipBranchSitePath(const Link* from);
	void SlaveWriteMappings();

	void GlobalSetTestData();
	void SlaveSetTestData();

	// private:

	int** siteprofilesuffstatcount;
	int* allocsiteprofilesuffstatcount;

	ZippedSequenceAlignment* zipdata;
	SequenceAlignment* truedata;
};

#endif
