
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


class MCParameters;
class Bipartition;
class PBTree;

class TaxaParameters	{

	public:

	TaxaParameters();
	TaxaParameters(string filename);
	TaxaParameters(int N, string* inNames = 0);
	TaxaParameters(MCParameters* inParam);
	TaxaParameters(PBTree* inTree);
	TaxaParameters(const TaxaParameters& from);
	~TaxaParameters();

	void				ReadFromStream(istream& is);
	void				WriteToStream(ostream& os);
	void				ReadNexus(istream& is);

	
	string	GetSpeciesName(int index);
	Boolean Member(int index);

	int		RegisterWith(TaxaParameters* inParam);
	MCParameters*	GetParameters()	{return mParam;}

	Bipartition			GetOutGroup();
	void				SetOutGroup(Bipartition& inOutGroup);

	MCParameters*			mParam;
	int				Ntaxa;				// internal number of taxa
	string*				SpeciesNames;
		
}
;
