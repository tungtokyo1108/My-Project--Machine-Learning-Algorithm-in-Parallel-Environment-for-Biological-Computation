
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


#include <iostream>
#include <cstdlib>
using namespace std;

#include "CodonSequenceAlignment.h"
// #include "Exception.h"

CodonSequenceAlignment::CodonSequenceAlignment(SequenceAlignment* from, bool force_stops, GeneticCodeType type)	{

		try	{
			DNAsource = from;

			if (from->Nsite % 3 != 0)	{
				cerr << "not multiple of three\n";
				exit(1);
			}		 
			Nsite = from->Nsite/3;
			Ntaxa = from->Ntaxa;
			statespace = new CodonStateSpace(type);

			taxset = DNAsource->GetTaxonSet();

			// make my own arrays
			// make translation
			Data = new int*[Ntaxa];
			for (int i=0; i<Ntaxa; i++)	{
				Data[i] = new int[Nsite];
				for (int j=0; j<Nsite; j++)	{
					try {
						Data[i][j] = GetCodonStateSpace()->GetCodonFromDNA(DNAsource->GetState(i, 3*j), DNAsource->GetState(i, 3*j+1), DNAsource->GetState(i, 3*j+2));
					}
					catch(...)	{
					// catch(Exception e)	{
						cerr << "in CodonSequenceAlignment: taxon " << i << " and codon " << j << " (site " << 3*j << ")\n";
						cerr << "taxon : " << taxset->GetTaxon(i) << '\n'; 
						if (force_stops)	{
							Data[i][j] = -1;
						}
						else	{
							throw;
						}
					}
				}
			}

		}
		catch(...)	{
		// catch(Exception)	{
			cerr << "Codon Sequence Alignment: failed to read the datafile\n";
			exit(1);
		}
	}


void CodonSequenceAlignment::ToStream(ostream& os)	{

	os << Ntaxa << '\t' << 3 * Nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}
	
	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			os << statespace->GetState(GetState(i,j));
		}
		os << '\n';
	}
	os << '\n';
}

/*void CodonSequenceAlignment::ToStream(ostream& os)	{

	os << Ntaxa << '\t' << Nsite << '\n';
	int max = 0;
	for (int i=0; i<Ntaxa; i++)	{
		int l = taxset->GetTaxon(i).length();
		if (max < l)	{
			max = l;
		}
	}
	
	for (int i=0; i<Ntaxa; i++)	{
		os << taxset->GetTaxon(i);
		for (unsigned int j=0; j< 5 + max - taxset->GetTaxon(i).length(); j++)	{
			os << ' ';
		}
		for (int j=0; j<Nsite; j++)	{
			int tempCodonState = Data[i][j];
			int tempAAState = GetCodonStateSpace()->Translation(tempCodonState);
			if (tempCodonState == unknown)	{
				os << '-';
			}
			else {
				os << AminoAcids[tempAAState];
			}
		}
		os << '\n';
	}
	os << '\n';
}*/
