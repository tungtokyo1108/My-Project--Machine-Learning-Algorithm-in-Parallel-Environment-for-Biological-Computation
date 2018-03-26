
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#ifndef ZIP_H
#define ZIP_H

#include "SequenceAlignment.h"

class ZippedSequenceAlignment : public SequenceAlignment	{

	public:

	ZippedSequenceAlignment(SequenceAlignment* infrom)	{

		from = infrom;

		Ntaxa = from->Ntaxa;
		Nsite = from->Nsite;
		taxset = from->taxset;
		statespace = from->statespace;
		Nstate = statespace->GetNstate();

		ComputeZipArrays();
	}

	~ZippedSequenceAlignment()	{
		// things here !!!
	}

	SequenceAlignment* GetTemplate() {return from;}
	int GetZipSize(int site) {return ZipSize[site];}
	int GetOrbitSize(int site) {return OrbitSize[site];}
	int GetStateFromZip(int site, int state) {
		if (state >= OrbitSize[site])	{
			cerr << "error in getstate from zip\n";
			exit(1);
		}
		return Indices[site][state];
	}
	bool InOrbit(int site, int state)	{
		return Orbit[site][state];
	}

	private:

	void ComputeZipArrays();

	SequenceAlignment* from;
	int Nstate;

	double SpeedFactor;

	bool** Orbit;
	int* OrbitSize;
	double MeanOrbitSize;
	int* ZipSize;
	int** Indices;
	int** ZipIndices;
	int** ZipData;

	
	

};


#endif

