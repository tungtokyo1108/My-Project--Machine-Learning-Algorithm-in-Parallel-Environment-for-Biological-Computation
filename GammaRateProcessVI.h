
/********************

**********************/


#ifndef GAMRATEVI_H
#define GAMRATEVI_H

#include "RateProcessVI.h"
// #include "BranchProcessVI.h"
#include <gsl/gsl_sf_psi.h>

class GammaRateProcessVI : public virtual RateProcessVI {

	public:

	GammaRateProcessVI() : Ncat(0), rate(0), RateAlphaVI(0), RateBetaVI(0), RateAlphahat(0), RateBetahat(0) {}
	virtual ~GammaRateProcessVI() {}
      
        /*void TotalParameterRate() {
             
             TotalParaRate = 0;
               for(int i=0; i<GetNsite(); i++) {
               TotalParaRate += (RateAlphaVI[i]/RateBetaVI[i]);
               }  
        }

        double TotalParaRate;
        double GetTotalParaRate() {return TotalParaRate;}

        //--------------------------------------
        void TotalParameterLength() {
            
            TotalParaLength = 0;
            for(int i=1; i<GetNbranch(); i++) {
               TotalParaLength += (branchalphaVI[i] / branchbetaVI[i]);
            }
        }

        double TotalParaLength;
        double GetTotalParaLength() {return TotalParaLength;}*/
        //--------------------------------------

	/*double GetAlpha() {
     
              double total = 0;
              for(int i=0; i<GetNsite(); i++) {
              // for(int i=GetSiteMin(); i<GetSiteMax(); i++) {
                     total += (GetRateAlphaVI(i)/GetRateBetaVI(i));
              }
              // return total/(GetSiteMax() - GetSiteMin());
              return total / GetNsite();
        }*/

        double GetAlpha()  {

               double total = 0;
               double* ratealphavi = new double[GetNsite()];
               double* ratebetavi = new double[GetNsite()];
               for (int i=0; i<GetNsite(); i++)  {
                       ratealphavi[i] = 0.0;
                       ratebetavi[i] = 0.0;
               } 
               GlobalGetRateAlphaVI();
               GlobalGetRateBetaVI();
               for (int i=0; i<GetNsite(); i++)  {
                       ratealphavi[i] += RateAlphaVI[i];
                       ratebetavi[i] += RateBetaVI[i];
                       total += ratealphavi[i]/ratebetavi[i]; 
               }
               return total;
        }
  
        /*double GetAlpha()   {
               
               double total = GlobalGetMeanRateVI();
               return total;
        }*/

	int GetNrate(int site)	{
		if (SumOverRateAllocations())	{
			return Ncat;
		}
		return 1;
	}

	int GetNcat() {return Ncat;}

	double GetRate(int site)	{
		// cat should be == 0
		/*if (SumOverRateAllocations())	{
			return rate[cat];
		}*/
		return rate[site];
	}

        virtual double GetRateAlphaVI(int site)  {
               return RateAlphaVI[site];
        }

        virtual double GetRateBetaVI(int site)  {
               return RateBetaVI[site];
        } 

	double GetRateWeight(int site, int cat)	{
		if (SumOverRateAllocations())	{
			return 1.0/Ncat;
		}
		return 1.0;
	}

	void ActivateSumOverRateAllocations() {
		sumflag = true;
	}

	void InactivateSumOverRateAllocations(int* ratealloc) {
		for (int i=0; i<GetNsite(); i++)	{
			alloc[i] = ratealloc[i];
		}
		sumflag = false;
	}

	double GetPriorMeanRate()	{
		double total = 0;
		for (int i=0; i<GetNsite(); i++)	{
			total += rate[i];
		}
		return total / GetNsite();
	}

        void SetAlpha(double inalpha)  {
             alpha = inalpha;
	     UpdateDiscreteCategories();
        }        

	void ToStream(ostream& os);
	void FromStream(istream& is);

	// protected:

	void Create(int innsite, int inncat);
	void Delete();

	void SampleRate();


	// void GlobalUpdateRateSuffStat();
	// void SlaveUpdateRateSuffStat();
	void UpdateRateSuffStat();
     
	double UpdateDiscreteCategories();

	// void TotalParameterRate();
        // double GetTotalParaRate() {return TotalParaRate;}

        double ELBORate();
        double ELBO_Rate;
        double GetRateAlphahat(int site) {return RateAlphahat[site];}
        double GetRateBetahat(int site)  {return RateBetahat[site];}
	int Ncat;
	double* rate;
	double beta;
        double alpha;
        double taurate;
        double kapparate;
        double* LearningRateRate;
        double* RateAlphahat;
        double* RateBetahat;
        double* RateAlphaVI;
        double* RateBetaVI; 
	int* alloc;
	int* ratesuffstatcount;
	double* ratesuffstatbeta;
        
        double LogRatePrior();

        void EstimateRateParameter();
        void SlaveEstimateRateParameter();
        void GlobalEstimateRateParameter();

        double MoveRateParameter();
        void SlaveMoveRateParameter();
        void GlobalMoveRateParameter();

        double MoveRate();
        void GlobalMoveRate();
        void SlaveMoveRate();

        double Move();

        void UpdateRateProcess();
        virtual double GlobalGetMeanRateVI();
        virtual void SlaveGetMeanRateVI();
        // double MeanRateVI;
        void GlobalGetRateAlphaVI();
        void SlaveGetRateAlphaVI();
        void GlobalGetRateBetaVI();
        void SlaveGetRateBetaVI();
        double GlobalELBORate();
        void SlaveELBORate();

};

#endif

