
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#ifndef POISSONPROFILE_H
#define POISSONPROFILE_H

#include "ProfileProcess.h"

// Poisson (F81) implementation
class PoissonProfileProcess : public virtual ProfileProcess {

	public:

	PoissonProfileProcess() {}
	virtual ~PoissonProfileProcess() {}

	protected:

	virtual void UpdateZip(int site) = 0;

	// implemented in specialized phyloprocess subclasses
	virtual const int* GetSiteProfileSuffStatCount(int site) = 0;

};

#endif

