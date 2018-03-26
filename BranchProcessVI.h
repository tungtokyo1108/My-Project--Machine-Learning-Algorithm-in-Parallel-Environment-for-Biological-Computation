
/********************

**********************/


#ifndef BRANCHVI_H
#define BRANCHVI_H

#include "Tree.h"
#include "Chrono.h"

class BranchProcessVI : public NewickTree {

	public:

	BranchProcessVI() : tree(0), blarray(0) {}
	virtual ~BranchProcessVI() {}

	virtual string GetVersion() = 0;

	// 
	NewickTree* GetLengthTree() {return this;}
	Tree* GetTree() {return tree;}
	Link* GetRoot() {return tree->GetRoot();}
	const Link* GetRoot() const {return tree->GetRoot();}

	string GetBranchName(const Link* link) const;
	string GetNodeName(const Link* link) const;

	double GetLength(const Branch* branch) {return branch ? blarray[branch->GetIndex()] : 0;}
	double GetLength(const Branch* branch) const {return branch ? blarray[branch->GetIndex()] : 0;}
        double GetBranchAlpha(const Branch* branch) {return branch ? branchalphaVI[branch->GetIndex()] : 0;}
        double GetBranchBeta(const Branch* branch) {return branch ? branchbetaVI[branch->GetIndex()] : 0;}

	void SetLength(const Branch* branch, double inlength) {
		if (! branch)	{
			cerr << "error in branch process: null branch\n";
			exit(1);
		}
		if (! blarray)	{
			cerr << "array not created\n";
			exit(1);
		}
		blarray[branch->GetIndex()] = inlength;
	}

	void Backup();
	void Restore();

	double ProposeMove(const Branch* branch, double tuning);
	void MoveBranch(const Branch* branch, double factor);
	void Restore(const Branch* branch);

	virtual double GetNormFactor()  = 0;

	// return number of branches
	int MoveAllBranches(double factor);
	int RecursiveMoveAllBranches(const Link* from, double e);

	virtual void SampleLength(const Branch* branch) = 0;
	virtual double LogBranchLengthPrior(const Branch* branch) = 0;
        virtual double GetbranchalphaVI(int index) = 0;
        virtual double GetbranchbetaVI(int index) = 0;

	double GetTotalLength()	{
		return RecursiveTotalLength(GetRoot());
	}

        double GetTotalLengthVI()	{
		return RecursiveTotalLengthVI(GetRoot());
	}
   
        /*double GetTotalLengthVI()  {
                return RecursiveTotalLengthVI();
        }*/

	double GetRenormTotalLength()	{
		// return RecursiveTotalLength(GetRoot());
                return RecursiveTotalLength(GetRoot()) * GetNormFactor(); 
	}

	void RenormalizeBranchLengths()	{
		double tmp = GetNormFactor();
                // double tmp = 1;
		RecursiveNormalizeBranchLengths(GetRoot(),tmp);
	}

	void DenormalizeBranchLengths()	{
		double tmp = GetNormFactor();
                // double tmp = 1;
		RecursiveNormalizeBranchLengths(GetRoot(),1.0 / tmp);
	}

	void RecursiveNormalizeBranchLengths(const Link* from, double factor);
	
	double RecursiveTotalLength(const Link* from);
        double RecursiveTotalLengthVI(const Link* from);
        // double RecursiveTotalLengthVI();

	void SetLengthsFromNames()	{
		if (blarray)	{
			RecursiveSetLengthsFromNames(GetRoot());
		}
		else	{
			cerr << "set lengths from names called without blarray\n";
		}
	}
	
	void SetNamesFromLengths()	{
		if (blarray)	{
			RecursiveSetNamesFromLengths(GetRoot());
		}
		else	{
			cerr << "set names from length called without blarray\n";
		}
	}

	void RecursiveSetLengthsFromNames(const Link* from);
	void RecursiveSetNamesFromLengths(const Link* from);

	void SampleLength()	{
		RecursiveSampleLength(GetRoot());
	}

	double LogLengthPrior()	{
		return RecursiveLogLengthPrior(GetRoot());
	}

	virtual void GlobalUpdateBranchLengthSuffStat() = 0;
	virtual void SlaveUpdateBranchLengthSuffStat() = 0;
	virtual void UpdateBranchLengthSuffStat() = 0;
        
        /*virtual void TotalParameterRate() = 0;
        virtual double GetTotalParaRate() = 0;*/

	/*
	virtual double GetBranchLengthSuffStatBeta(const Branch* branch) = 0;
	virtual int GetBranchLengthSuffStatCount(const Branch* branch) = 0;
	*/

	virtual double GetBranchLengthSuffStatBeta(int index) = 0;
        // virtual int GetBranchLengthSuffStatCountELBO(int index) = 0;
	virtual int GetBranchLengthSuffStatCount(int index) = 0;

	// implements a map<const Branch*, double>

	// Move function ?
	// how about the tuning parameters ?

	double LengthSuffStatLogProb();

	virtual void ToStream(ostream& os) = 0;
	virtual void FromStream(istream& is) = 0;

	virtual const TaxonSet* GetTaxonSet() const = 0;
        virtual int GetMyid() = 0;
        virtual int GetNprocs() = 0;

	// protected:

	int GetNbranch()	{
		return tree->GetNbranch();
	}

	int GetNnode()	{
		return tree->GetNnode();
	}

	int GetNlink()	{
		return tree->GetNlink();
	}

	virtual void Create(Tree* intree) {
		tree = intree;
		blarray = new double[GetNbranch()];
		bkarray = new double[GetNbranch()];
		blarray[0] = 0;
                LearningRateBranch = new double[GetNbranch()];
                branchalphaVI = new double[GetNbranch()];
                branchbetaVI = new double[GetNbranch()];
                branchalphahat = new double[GetNbranch()];
                branchbetahat = new double[GetNbranch()];
                for (int i=1; i<GetNbranch(); i++)	{
                     // branchalphaVI[i] = rnd::GetRandom().sExpo();
                     // branchbetaVI[i] = rnd::GetRandom().sExpo();
                     branchalphaVI[i] = 1;
                     branchbetaVI[i] = 10;
                }
		SetLengthsFromNames();
	}
	virtual void Delete() {
		delete[] blarray;
		delete[] bkarray;
                delete[] branchalphaVI;
                delete[] branchbetaVI;
                delete[] LearningRateBranch; 
	}

	double RecursiveLogLengthPrior(const Link* from);
	void RecursiveSampleLength(const Link* from);

	Tree* tree;
	double* blarray;
	double* bkarray;

        double* branchalphaVI;
	double* branchbetaVI;
        double* branchalphahat;
        double* branchbetahat;
        double* LearningRateBranch;

        /*virtual void TotalParameterLength() = 0;

        double TotalParaLength;
        virtual double GetTotalParaLength() = 0;       
        */

	Chrono chronolength;

};

#endif

