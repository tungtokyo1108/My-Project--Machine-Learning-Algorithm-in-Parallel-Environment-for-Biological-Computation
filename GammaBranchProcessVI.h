
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#ifndef GAMMABRANCHVI_H
#define GAMMABRANCHVI_H

#include "TaxonSet.h"
#include "BranchProcessVI.h"
#include "RateProcessVI.h"
#include "Random.h"

class GammaBranchProcessVI : public virtual BranchProcessVI {

	public:

	GammaBranchProcessVI() : betaprior(0) {}
	virtual ~GammaBranchProcessVI() {}
        
	double LogBranchLengthPrior(const Branch* branch);

        double DQBranchAlpha(const Branch* branch);
        double DQBranchBeta(const Branch* branch);   
        double ELBOLength();
        double ELBO_Length;
        double GlobalELBOLength();
        void SlaveELBOLength();       
        void EstimateBranchParameter();
        void SlaveEstimateBranchParameter();
        void GlobalEstimateBranchParameter();

        double MoveBranchParamter();
        void SlaveMoveBranchParamter();
        void GlobalMoveBranchParamter();

        double MoveLength();
        void SlaveMoveLength();
        void GlobalMoveLength();

        double Move();

	void SampleLength();
	void SampleLength(const Branch* branch);
       
	void ToStreamWithLengths(ostream& os, const Link* from);

	void ToStream(ostream& os);
	void FromStream(istream& is);

	protected:

	virtual void Create(Tree* intree, double inalpha = 1, double inbeta = 10)	{
		BranchProcessVI::Create(intree);
		// branchalpha = rnd::GetRandom().sExpo();
		// branchbeta = rnd::GetRandom().sExpo();
                branchalpha = 1;
                branchbeta = 10;
                taubranch = 1;
                kappabranch = 0.5;
	}

	void Delete() {}
        void TotalParameterLength();
        double GetTotalParaLength() {return TotalParaLength;}
        double TotalParaLength;        
        double branchalpha;
        double branchbeta;
        double taubranch;
        double kappabranch;
        double Getbranchalphahat(int index) {return branchalphahat[index];}
        double Getbranchbetahat(int index)  {return branchbetahat[index];}
        virtual double GetbranchalphaVI(int index)  {return branchalphaVI[index];}
        virtual double GetbranchbetaVI(int index)  {return branchbetaVI[index];}
	int betaprior;
        double dqbranchalpha;
        
};

#endif

