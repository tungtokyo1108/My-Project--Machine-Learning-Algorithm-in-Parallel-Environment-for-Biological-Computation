
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


class TaxaParameters;
class PBTree;
class PolyNode;
class Bipartition;
class PhyloBayes;

class TreeList	{


	public:

				TreeList(TaxaParameters* inParam, int inSize);
				TreeList();
				TreeList(string filename);
				~TreeList();

	PBTree*			GetTree(int index);

	int			GetSize()	{return mSize;}

	TaxaParameters*		GetParameters()	const {return mParam;}
	void			SetParameters(); // set parameters for each tree of the list 

	void			RootAt(Bipartition outgroup);

	void			ReadFromStream(istream& is);
	void			WriteToStream(ostream& os, int header = 0, int withLengths = 1, int withProbs = 0, int withSpeciesNames = 1,int withInternalLabels = 0);
	void			ToPS(string target, int every = 1, double sizeX=12, double sizeY=20, int withLengths=0, int withProbs = 0, int withSpeciesNames = 1, int withInternalLabels = 0);

	TaxaParameters*		mParam;
	int 			mSize;
	PBTree**		mTreeArray;

}
;


