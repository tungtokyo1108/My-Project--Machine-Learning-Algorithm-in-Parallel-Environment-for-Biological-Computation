
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/

class PolyNode;
class BipartitionList;
class Bipartition;
class TreeList;

class Consensus		:   public PBTree	{

	public:

	Consensus(TreeList* inTreeList, double* probarray = 0, double cutoff = 0.5);
	Consensus(BipartitionList* inBList, double cutoff = 0.5);

	~Consensus();

	BipartitionList*				MakeConsensus(double cutoff);

	BipartitionList*				GetBList()		const   {return mBList;}
	TreeList*					GetTreeList()	const {return mTreeList;}
	
	void						Insert( Bipartition& inPartition, double prob, double length);
	TreeList*					mTreeList;
	BipartitionList*				mBList;
	BipartitionList*				mReducedList;

}
;
