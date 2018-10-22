
/********************

**********************/


#include "BranchProcessVI.h"
#include "Random.h"

#include <sstream>
#include <iomanip>

string BranchProcessVI::GetNodeName(const Link* link)  const	{
	if (link->isLeaf())	{
		return link->GetNode()->GetName();
	}
	else	{
		return "";
	}
}

string BranchProcessVI::GetBranchName(const Link* link) const	{
	if (link->isRoot())	{
		return "";
	}
	else	{
		ostringstream s;
		s << GetLength(link->GetBranch());
		return s.str();
	}
}	

void BranchProcessVI::Backup()	{
	for (int i=0; i<GetNbranch(); i++)	{
		bkarray[i] = blarray[i];
	}
}

void BranchProcessVI::Restore()	{
	for (int i=0; i<GetNbranch(); i++)	{
		blarray[i] = bkarray[i];
	}
}

void BranchProcessVI::Restore(const Branch* branch)	{
	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	blarray[index] = bkarray[index];
}

void BranchProcessVI::MoveBranch(const Branch* branch, double m)	{

	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	bkarray[index] = blarray[index];
	blarray[index] *= exp(m);
}

double BranchProcessVI::ProposeMove(const Branch* branch, double tuning)	{

	if (! branch)	{
		cerr << "error in branchprocess::Movebranch: null branch\n";
		exit(1);
	}
	int index = branch->GetIndex();
	bkarray[index] = blarray[index];
	double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
	blarray[index] *= exp(m);
	return m;
}

int BranchProcessVI::MoveAllBranches(double e)	{
	return RecursiveMoveAllBranches(GetRoot(),e);
}

int BranchProcessVI::RecursiveMoveAllBranches(const Link* from, double e)	{

	int n = 0;
	if (!from->isRoot())	{
		int index = from->GetBranch()->GetIndex();
		blarray[index] *= e;
		n++;
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		n += RecursiveMoveAllBranches(link->Out(),e);
	}
	return n;
}

double BranchProcessVI::RecursiveLogLengthPrior(const Link* from)	{

	double total = 0;
	if (! from->isRoot())	{
		total += LogBranchLengthPrior(from->GetBranch());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveLogLengthPrior(link->Out());
	}
	return total;
}

void BranchProcessVI::RecursiveSampleLength(const Link* from)	{

	if (! from->isRoot())	{
		SampleLength(from->GetBranch());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSampleLength(link->Out());
	}
}

void BranchProcessVI::RecursiveNormalizeBranchLengths(const Link* from, double factor)	{

	if (! from->isRoot())	{
		SetLength(from->GetBranch(),GetLength(from->GetBranch())*factor);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveNormalizeBranchLengths(link->Out(),factor);
	}
}

double BranchProcessVI::RecursiveTotalLength(const Link* from)	{
	double total = 0;
	if (! from->isRoot())	{
		total += GetLength(from->GetBranch());
	}
	else	{
		if (GetLength(from->GetBranch()))	{
			cerr << "error: non null branch length for root\n";
			cerr << GetLength(from->GetBranch()) << '\n';
			exit(1);
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += RecursiveTotalLength(link->Out());
	}
	return total;
}

double BranchProcessVI::RecursiveTotalLengthVI(const Link* from) {

        double total = 0;
        if (!from->isRoot()) {
                 int index = from->GetBranch()->GetIndex();
                 total += branchalphaVI[index] / branchbetaVI[index];
        }
        for (const Link* link=from->Next(); link!=from; link=link->Next()) {
                  total += RecursiveTotalLengthVI(link->Out());
        }
        return total;
}

void BranchProcessVI::RecursiveSetLengthsFromNames(const Link* from)	{
	if (! from->isRoot())	{
		double l = atof(from->GetBranch()->GetName().c_str());
		if (l <= 0)	{
			l = 0.1;
			/*
			cerr << "error in BranchProcess::SetLengthsFromFile : " << l << '\n';
			exit(1);
			*/
		}
		SetLength(from->GetBranch(),l);
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetLengthsFromNames(link->Out());
	}
}

void BranchProcessVI::RecursiveSetNamesFromLengths(const Link* from)	{

	if (! from->isRoot())	{
		ostringstream s;
		s.precision(15);
		s << GetLength(from->GetBranch());
		from->GetBranch()->SetName(s.str());
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetNamesFromLengths(link->Out());
	}
}

double BranchProcessVI::LengthSuffStatLogProb()	{

	double total = 0;
	for (int i=1; i<GetNbranch(); i++)	{
		total += GetBranchLengthSuffStatCount(i) * log(blarray[i]);
		total -= GetBranchLengthSuffStatBeta(i) * blarray[i];
	}
	return total;
}



