
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#ifndef MPIMODULE_H
#define MPIMODULE_H

class MPIModule {

	protected:

	MPIModule() {}
	virtual ~MPIModule() {}

	virtual int GetNprocs() = 0;
	virtual int GetMyid() = 0;
	virtual int GetSiteMin(int proc = -1) = 0;
	virtual int GetSiteMax(int proc = -1) = 0;
	virtual void MakeMPIPartition() = 0;

};


#endif

