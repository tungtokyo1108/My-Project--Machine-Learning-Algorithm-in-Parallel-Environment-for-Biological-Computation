
/********************


**********************/


#ifndef RASCATVI_H
#define RASCATVI_H

#include "PoissonSubstitutionProcess.h"
#include "PoissonPhyloProcess.h"
#include "GammaRateProcessVI.h"
#include "PoissonDPProfileProcessVI.h"
#include "GammaBranchProcessVI.h"
#include "PoissonSBDPProfileProcessVI.h"

class RASCATSubstitutionProcessVI : public virtual PoissonSubstitutionProcess, public virtual GammaRateProcessVI, public virtual PoissonDPProfileProcessVI, public virtual PoissonSBDPProfileProcessVI {

	public:

	RASCATSubstitutionProcessVI() {}
	virtual ~RASCATSubstitutionProcessVI() {}

	protected:

	virtual void Create(int, int)	{
		cerr << "error : in RASCATSubProcess::Create(int,int)\n";
		exit(1);
	}

	virtual void Create(int nsite, int ncat, int nstate,int insitemin,int insitemax)	{
		PoissonSubstitutionProcess::Create(nsite,nstate,insitemin,insitemax);
		GammaRateProcessVI::Create(nsite,ncat);
		PoissonDPProfileProcessVI::Create(nsite,nstate);
	}

	virtual void Delete()	{
		PoissonDPProfileProcessVI::Delete();
		GammaRateProcessVI::Delete();
		PoissonSubstitutionProcess::Delete();
	}

};

class RASCATGammaPhyloProcessVI : public virtual PoissonPhyloProcess, public virtual RASCATSubstitutionProcessVI, public virtual GammaBranchProcessVI	{

	public:

        virtual void SlaveExecute(MESSAGE);
	virtual void GlobalUpdateParameters();
	virtual void SlaveUpdateParameters();

	RASCATGammaPhyloProcessVI() {}

	RASCATGammaPhyloProcessVI(string indatafile, string treefile, int nratecat, int iniscodon, GeneticCodeType incodetype, int infixtopo, int inkappaprior, double inmintotweight, int indc, int me, int np)	{
		myid = me;
		nprocs = np;

		fixtopo = infixtopo;
		dc = indc;
		kappaprior = inkappaprior;
		SetMinTotWeight(inmintotweight);

		datafile = indatafile;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
		if (dc)	{
			plaindata->DeleteConstantSites();
		}
		const TaxonSet* taxonset = plaindata->GetTaxonSet();
		if (treefile == "None")	{
			tree = new Tree(taxonset);
			if (myid == 0)	{
				tree->MakeRandomTree();
				GlobalBroadcastTree();
			}
			else	{
				SlaveBroadcastTree();
			}
		}
		else	{
			tree = new Tree(treefile);
		}
		tree->RegisterWith(taxonset,myid);
		
		int insitemin = -1,insitemax = -1;
		if (myid > 0) {
			int width = plaindata->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = plaindata->GetNsite();
			}
			else {
				insitemax = myid*width;
			}
		}

		Create(tree,plaindata,nratecat,insitemin,insitemax);

		if (myid == 0)	{
			Sample();
			GlobalUnfold();
		}
	}

	RASCATGammaPhyloProcessVI(istream& is, int me, int np)	{
		myid = me;
		nprocs = np;

		FromStreamHeader(is);
		is >> datafile;
		int nratecat;
		is >> nratecat;
		if (atof(version.substr(0,3).c_str()) > 1.3)	{
			is >> iscodon;
			is >> codetype;
			is >> kappaprior;
			is >> mintotweight;
		}
		else	{
			iscodon = 0;
			codetype = Universal;
			kappaprior = 0;
			mintotweight = -1;
		}
		is >> fixtopo;
		if (atof(version.substr(0,3).c_str()) > 1.4)	{
			is >> NSPR;
			is >> NNNI;
		}
		else	{
			NSPR = 10;
			NNNI = 0;
		}
		is >> dc;
		SequenceAlignment* plaindata = new FileSequenceAlignment(datafile,0,myid);
		if (dc)	{
			plaindata->DeleteConstantSites();
		}
		const TaxonSet* taxonset = plaindata->GetTaxonSet();

		int insitemin = -1,insitemax = -1;
		if (myid > 0) {
			int width = plaindata->GetNsite()/(nprocs-1);
			insitemin = (myid-1)*width;
			if (myid == (nprocs-1)) {
				insitemax = plaindata->GetNsite();
			}
			else {
				insitemax = myid*width;
			}
		}

		tree = new Tree(taxonset);
		if (myid == 0)	{
			tree->ReadFromStream(is);
			GlobalBroadcastTree();
		}
		else	{
			SlaveBroadcastTree();
		}
		tree->RegisterWith(taxonset,0);

		Create(tree,plaindata,nratecat,insitemin,insitemax);

		if (myid == 0)	{
			FromStream(is);
			GlobalUnfold();
		}
	}

	~RASCATGammaPhyloProcessVI() {
		Delete();
	}

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		// yet to be implemented
		return 0;
	}

	double GetLogLikelihood()	{
		return logL;
	}

        double ELBO;
        double pre_lower_bound;
        double lower_bound;
        double change;
        double GetELBO () {
              return ELBO;
        }
        // double CheckELBO();        

	void TraceHeader(ostream& os)	{
		os << "#iter\ttime\ttopo\tloglik\tlength\talphaVI\tNmode\tstatent";
                // os << "\ELBO";
		os << "\tkappaVI\tallocent";
		os << '\n'; 
	}

	void Trace(ostream& os)	{

		os << GetSize();
		if (chronototal.GetTime())	{
			os << '\t' << chronototal.GetTime() / 1000;
			os << '\t' << ((int) (propchrono.GetTime() / chronototal.GetTime() * 100));
			chronototal.Reset();
			propchrono.Reset();
		}
		else	{
			os << '\t' << 0;
			os << '\t' << 0;
		}

		os << '\t' << GetLogLikelihood() << '\t' << GetRenormTotalLength() << '\t' << GetAlpha();
		os << '\t' << GetNOccupiedComponent() << '\t' << GetStatEnt();
                os << '\t' << change;
		// os << '\t' << GetMeanDirWeight();
		// os << '\t' << kappa_alpha << '\t' << GetAllocEntropy();

		os << '\n';
	}

	virtual double Move(double tuning = 1.0)	{

		chronototal.Start();
		propchrono.Start();
		BranchLengthMove(tuning);
		BranchLengthMove(0.1 * tuning);
		if (! fixtopo)	{
			MoveTopo(10,0);
		}
		propchrono.Stop();

		GlobalCollapse();

		GammaBranchProcessVI::Move();

		// this one is important 
		GlobalUpdateParameters();
		GammaRateProcessVI::Move();
		
		// RASCATSubstitutionProcess::MoveRate(tuning);

		// this one is not useful
		// because uniformized process:
		// conditional on discrete substitution mapping
		// profiles do not depend on branch lengths and site rates
		// GlobalUpdateParameters();

		PoissonDPProfileProcessVI::Move(1,1,1);

		GlobalUnfold();
		chronototal.Stop();

		// Trace(cerr);

		return 1;
	
	}

	virtual void ReadPB(int argc, char* argv[]);
	void ReadSiteProfiles(string name, int burnin, int every, int until);

	void ToStreamHeader(ostream& os)	{
		PhyloProcess::ToStreamHeader(os);
		os << datafile << '\n';
		os << GetNcat() << '\n';
		os << iscodon << '\n';
		os << codetype << '\n';
		os << kappaprior << '\n';
		os << mintotweight << '\n';
		os << fixtopo << '\n';
		os << NSPR << '\t' << NNNI << '\n';
		os << dc << '\n';
		SetNamesFromLengths();
		GetTree()->ToStream(os);
	}

	void ToStream(ostream& os)	{
		GammaBranchProcessVI::ToStream(os);
		GammaRateProcessVI::ToStream(os);
		PoissonDPProfileProcessVI::ToStream(os);
	}

	void FromStream(istream& is)	{
		GammaBranchProcessVI::FromStream(is);
		GammaRateProcessVI::FromStream(is);
		PoissonDPProfileProcessVI::FromStream(is);
		GlobalUpdateParameters();
	}


	virtual void Create(Tree* intree, SequenceAlignment* indata, int ncat,int insitemin,int insitemax)	{
		PoissonPhyloProcess::Create(intree,indata);
		// PoissonPhyloProcess::Create(intree,indata,indata->GetNstate(),insitemin,insitemax);
		RASCATSubstitutionProcessVI::Create(indata->GetNsite(),ncat,indata->GetNstate(),insitemin,insitemax);
		GammaBranchProcessVI::Create(intree);
	}
		
	virtual void Delete()	{
		GammaBranchProcessVI::Delete();
		RASCATSubstitutionProcessVI::Delete();
		PoissonPhyloProcess::Delete();
	}

	int iscodon;
	GeneticCodeType codetype;
	int fixtopo;
	int NSPR;
	int NNNI;
	int dc;
};

#endif

