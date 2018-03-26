
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/


class Bipartition;
class TreeList;
class Consensus;
class PBTree;

class PolyNode	{

	public:

	PolyNode();
	PolyNode(int inLabel, int inLevel = -1);

	~PolyNode();
	// deletion is not recursive

	void		deleteRecursive();
	// accessors

	PolyNode*	Up()	{return up;}
	PolyNode*	Down()	{return down;}
	PolyNode*	Next()	{return next;}
	PolyNode*	Prev()	{return prev;}

	int			IsLeaf() {return (down == 0);}
	int			IsRoot() {return (up == 0);}

	int			GetDegree() {return degree;}
	int			GetMaxDegree();
	void			ComputeDegree();
	int			CheckDegree();
	int			IsDichotomous();
	void			GetLeafSet(Bipartition& bp);

	int 			GetLabel() {return label;}
	void			SetLabel(int inLabel)	{ label = inLabel;}

	string			GetName()	{return name;}
	void			SetName(string inName)	{ name = inName;}

	void			SetNames();

	double			GetBranchLength()	{return branchLength;}
	void			SetBranchLength(double inBranchLength)	{branchLength = inBranchLength;}

	double			GetProb()	{return mProb;}
	void			SetProb(double inProb)	{mProb = inProb;}

	void			SetTree(PBTree* inTree)	{mTree = inTree;}
	PBTree*			GetTree()	{return mTree;}

	PolyNode*		Simplify();
	void			Detach();
	// void			DetachAndRemove(); // this one also removes the up node, if it ends up with only one descendant
	void			DetachRecursive();
	PolyNode*		AttachTo(PolyNode* inNode);

	int 			ComputeMinLeafLabel();

	void			DetachOld(PolyNode* inNode);
	void			Insert(PolyNode* inNode);	// will be inserted as the last one


	int			RegisterWithData(int& currentLabel, int* found, string* SpeciesNames, int Ntaxa);
	int			Dichotomise();
	void			SetSuper(PBTree* inTree);			// recursive 

	int			GetSize();
	double			GetDepth();
	double			GetMinDepth();
	double			GetLength();
	int			GetIntDepth();
	void			ComputeSizeAndDepth();
	void			SetFormalBranchLengths(double scale, int order);

	void			GetSpeciesNames(string* name, int& index);
	void			GetSpeciesNames(string* name);

	PolyNode*		FindNode(Bipartition& inPartition, Boolean& Orientation);

	PolyNode*		FindNodePruning(Bipartition& inPartition,
										Boolean& Orientation,
										Bipartition& outPartition);

	Bipartition		GetBipartition();
	Bipartition		BipartitionPruning(BipartitionList* inBList);
	Bipartition		BipartitionPruningWithSupports(BipartitionList* inBList);

	int			Analyse(Bipartition& inPartition);
	void			Insert(Bipartition& inPartition, double prob, double length);

	PolyNode*		FlipFlop();

	void			TranslateLabels(int offset);
	int			SortLeaves();
	string			SortLeavesAlphabetical();
	void 			Phylip(ostream& os, int withLengths = 0, int withProbs = 1, int withSpeciesNames = 1, int withInternalLabels = 0);
	void			ResetLabels();
	/*
	void			GetRootedBPList(int** bplist, int Nnode, int Ntaxa);
	void			Revert();
	void			ReverseMapLabels(int* reversemap);
	*/

	PBTree* 	mTree;

	PolyNode* up;
	PolyNode* down;
	PolyNode* next;
	PolyNode* prev;


	int label;
	int degree;

	int value;

	string name;

	double branchLength;
	double mProb;

	int 	size;
	int	intdepth;
	double	depth;

	double X;
	double Y;
	double Y2;
	void Swap();
}
;
