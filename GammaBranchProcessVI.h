
/********************


**********************/


#ifndef GAMMABRANCHVI_H
#define GAMMABRANCHVI_H

#include "TaxonSet.h"
#include "BranchProcessVI.h"
#include "RateProcessVI.h"
#include <gsl/gsl_sf_psi.h>
#include "Random.h"

class GammaBranchProcessVI : public virtual BranchProcessVI {

	public:

	GammaBranchProcessVI() : betaprior(0) {}
	virtual ~GammaBranchProcessVI() {}
        
	double LogBranchLengthPrior(const Branch* branch);

        /*void TotalParameterLength() {
            
            TotalParaLength = 0;
            for(int i=1; i<GetNbranch(); i++) {
               TotalParaLength += (branchalphaVI[i] / branchbetaVI[i]);
            }
        }

        double TotalParaLength;
        double GetTotalParaLength() {return TotalParaLength;}

        //------------------------------------------------------
        void TotalParameterRate() {
             
             TotalParaRate = 0;
               for(int i=0; i<GetNsite(); i++) {
               TotalParaRate += (RateAlphaVI[i]/RateBetaVI[i]);
               }  
        }

        double TotalParaRate;
        double GetTotalParaRate() {return TotalParaRate;}*/
        //------------------------------------------------------

        double DQBranchAlpha(const Branch* branch);
        double DQBranchBeta(const Branch* branch);   
        double ELBOLength();
        double ELBO_Length;
        double GlobalELBOLength();
        void SlaveELBOLength();
     
        /*double GetbranchalphaVI(int GetNbranch()) {return branchalphaVI[GetNbranch()];}
        double GetbranchbetaVI(int GetNbramch())  {return branchbetaVI[GetNbranch()];}*/        

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
                /*LearningRateBranch = new double[GetNbranch()];
                branchalphaVI = new double[GetNbranch()];
                branchbetaVI = new double[GetNbranch()];
                branchalphahat = new double[GetNbranch()];
                branchbetahat = new double[GetNbranch()];
                for (int i=1; i<GetNbranch(); i++)	{
                     // branchalphaVI[i] = rnd::GetRandom().sExpo();
                     // branchbetaVI[i] = rnd::GetRandom().sExpo();
                     branchalphaVI[i] = 1;
                     branchbetaVI[i] = 10;
                }*/
		// RecursiveSampleLength(GetRoot());
	}

	void Delete() {}
        void TotalParameterLength();
        double GetTotalParaLength() {return TotalParaLength;}
        double TotalParaLength;        

	/*double* branchalphaVI;
	double* branchbetaVI;
        double* branchalphahat;
        double* branchbetahat;
        double* LearningRateBranch;*/
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

