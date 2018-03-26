
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#include "ProfileProcess.h"
#include "Random.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* ProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void ProfileProcess::Create(int innsite, int indim)	{
	if (nsite || dim)	{
		if (nsite != innsite)	{
			cerr << "error in phyloprocess creation: non matching number of sites\n";
			exit(1);
		}
		if (dim != indim)	{
			cerr << "error in phyloprocess creation: non matching number of states\n";
			cerr << dim << '\t' << indim << '\n';
			exit(1);
		}
	}
	else	{
		nsite = innsite;
		dim = indim;
	}
}

