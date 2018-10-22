
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#include "GammaRateProcessVI.h"
#include "Random.h"
#include "IncompleteGamma.h"

#include <cassert>
#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* GammaRateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void GammaRateProcessVI::Create(int innsite, int inncat)	{
	if (! rate)	{
		RateProcessVI::Create(innsite);
                Ncat = 4;
                // alpha = rnd::GetRandom().sExpo();
                // beta = rnd::GetRandom().sExpo();
                alpha = 1;
                beta = 1;
                taurate = 1;
                kapparate = 0.5;
                LearningRateRate = new double[GetNsite()];
                RateAlphahat = new double[GetNsite()];
                RateBetahat = new double[GetNsite()];
                RateAlphaVI = new double[GetNsite()];
                RateBetaVI = new double[GetNsite()];
                for(int i=0; i<GetNsite(); i++) {
                     // RateAlphaVI[i] = rnd::GetRandom().sExpo();
                     // RateBetaVI[i] = rnd::GetRandom().sExpo();
                     RateAlphaVI[i] = 1;
                     RateBetaVI[i] = 1;
                }
                rate = new double[GetNsite()];
		alloc = new int[GetNsite()];
		ratesuffstatcount = new int[GetNcat()];
		ratesuffstatbeta = new double[GetNcat()];
		// SampleRate();
	}
}


void GammaRateProcessVI::Delete() 	{
	delete[] rate;
	delete[] alloc;
	delete[] ratesuffstatcount;
	delete[] ratesuffstatbeta;
        delete[] RateAlphahat;
        delete[] RateBetahat;
        delete[] LearningRateRate; 
	rate = 0;
	alloc = 0;
	ratesuffstatcount = 0;
	ratesuffstatbeta = 0;
        RateProcessVI::Delete();
}

void GammaRateProcessVI::ToStream(ostream& os)	{
	os << alpha << '\n';
        
}

void GammaRateProcessVI::FromStream(istream& is)	{
        is >> alpha;
	SetAlpha(alpha);
}


// ------------------------------ Computate the variational distribution of GammaRateProcess -------------------------------------------------------------------

void GammaRateProcessVI::SampleRate()  {
     // alpha = rnd::GetRandom().sExpo();
     for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
     // for(int i=0; i<GetNsite(); i++) {
		rate[i] = rnd::GetRandom().Gamma(alpha, beta);
     }
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

double GammaRateProcessVI::UpdateDiscreteCategories()  {

        double* x = new double[GetNcat()];
	double* y = new double[GetNcat()];
	double lg = rnd::GetRandom().logGamma(alpha+1.0);
	for (int i=0; i<GetNcat(); i++)	{
		x[i] = PointGamma((i+1.0)/GetNcat(),alpha,alpha);
	}
	for (int i=0; i<GetNcat()-1; i++)	{
		y[i] = IncompleteGamma(alpha*x[i],alpha+1,lg);
	}
	y[GetNcat()-1] = 1.0;
	rate[0] = GetNcat() * y[0];
	for (int i=1; i<GetNcat(); i++)	{
		rate[i] = GetNcat() * (y[i] - y[i-1]);
	}
	delete[] x;
	delete[] y;
   return 1.0;   
}

double GammaRateProcessVI::LogRatePrior()	{
        
        double totalrate = 0; 
        for(int i=0; i<GetNsite(); i++) {
		totalrate += RateAlphaVI[i] * log(RateBetaVI[i]) - rnd::GetRandom().logGamma(RateAlphaVI[i]) + (RateAlphaVI[i]-1) * log(rate[i]) - RateBetaVI[i] * rate[i];
        }
	return totalrate;
}

double GammaRateProcessVI::ELBORate()   {

       ELBO_Rate = 0;
       // for (int i=0; i<GetNsite(); i++)  {
       for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
           ELBO_Rate += GetSiteRateSuffStatCount(i) * ((gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))));
           ELBO_Rate += (alpha * log(alpha) - rnd::GetRandom().logGamma(alpha) + (alpha - 1) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))) - alpha * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
           ELBO_Rate -= (GetRateAlphaVI(i) * log(GetRateBetaVI(i)) - rnd::GetRandom().logGamma(GetRateAlphaVI(i)));
           ELBO_Rate -= ((GetRateAlphaVI(i) - 1) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))) - GetRateBetaVI(i) * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
       }
       return ELBO_Rate;
}

// ---------------------------- Compute sufficient statistics of Rate Process ---------------------------------------------------------------------------------

/*void GammaRateProcessVI::GlobalUpdateRateSuffStat()	{
	assert(GetMyid() == 0);
	// MPI2
	// should ask the slaves to call their UpdateRateSuffStat
	// and then gather the statistics;
	int i,j,nprocs = GetNprocs(),workload = GetNcat();
	MPI_Status stat;
	MESSAGE signal = UPDATE_RATE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0; i<workload; ++i) {
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}
#ifdef BYTE_COM
	int k,l;
	double x;
	unsigned char* bvector = new unsigned char[workload*(sizeof(int)+sizeof(double))];

	for(i=1; i<nprocs; ++i) {
		MPI_Recv(bvector,workload*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<workload; ++j) {
			l = 0;
			for(k=sizeof(int)-1; k>=0; --k) {
				l = (l << 8) + bvector[sizeof(int)*j+k]; 
			}
			ratesuffstatcount[j] += l;
		}
		for(j=0; j<workload; ++j) {
			memcpy(&x,&bvector[sizeof(int)*workload+sizeof(double)*j],sizeof(double));
			ratesuffstatbeta[j] += x;
		}
	}
	delete[] bvector;
#else
	int ivector[workload];
	double dvector[workload];
        for(i=1; i<nprocs; ++i) {
                MPI_Recv(ivector,workload,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(j=0; j<workload; ++j) {
                        ratesuffstatcount[j] += ivector[j];                      
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(i=1; i<nprocs; ++i) {
                MPI_Recv(dvector,workload,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
                for(j=0; j<workload; ++j) {
                        ratesuffstatbeta[j] += dvector[j]; 
                }
        }
#endif
}
*/
void GammaRateProcessVI::UpdateRateSuffStat()	{

	for (int i=0; i<GetNcat(); i++)	{
		ratesuffstatcount[i] = 0;
		ratesuffstatbeta[i] = 0.0;
	}
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		ratesuffstatcount[alloc[i]] += GetSiteRateSuffStatCount(i);
		ratesuffstatbeta[alloc[i]] += GetSiteRateSuffStatBeta(i);
	}

}	

/*void GammaRateProcessVI::SlaveUpdateRateSuffStat()	{
	assert(GetMyid() > 0);

	UpdateRateSuffStat();

#ifdef BYTE_COM
	int n = 0;
	unsigned int j;
	unsigned char el_int[sizeof(int)],el_dbl[sizeof(double)];
	unsigned char* bvector = new unsigned char[GetNcat()*(sizeof(int)+sizeof(double))];

	for(int i=0; i<GetNcat(); ++i) {
		convert(el_int,ratesuffstatcount[i]);
		for(j=0; j<sizeof(int); ++j) {
			bvector[n] = el_int[j]; n++;
		}
	}
	for(int i=0; i<GetNcat(); ++i) {
		convert(el_dbl,ratesuffstatbeta[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_dbl[j]; n++;
		}
	}
	MPI_Send(bvector,GetNcat()*(sizeof(int)+sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;
#else
	MPI_Send(ratesuffstatcount,GetNcat(),MPI_INT,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(ratesuffstatbeta,GetNcat(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
#endif
}
*/

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                            Stochastic Optimization for Variational Parameters
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/*void GammaRateProcessVI::TotalParameterRate() {
             
             TotalParaRate = 0;
               for(int i=0; i<GetNsite(); i++) {
               TotalParaRate += (RateAlphaVI[i]/RateBetaVI[i]);
               }  
}*/

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                          Parallel Computation (Open-MPI) for Stochastic Optimization for Variational Parameters
// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double GammaRateProcessVI::Move() {
 
     GlobalUpdateSiteRateSuffStat();
     GlobalEstimateRateParameter();
     // EstimateRateParameter();
     GlobalMoveRateParameter();
     // MoveRateParameter();
     // MoveRate();
     GlobalMoveRate();
     // ELBORate();
     GlobalELBORate();
     return 1.0;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GammaRateProcessVI::SlaveGetRateAlphaVI()  {
 
        assert(GetMyid() > 0);
        
        int sitemin, sitemax, width = GetNsite()/(GetNprocs()-1);
        sitemin = (GetMyid()-1)*width;
        if (GetMyid() == (GetNprocs()-1))  {
                     sitemax = GetNsite();
        }
        else    {
                     sitemax = GetMyid()*width;
        }
        MPI_Send(RateAlphaVI+sitemin,sitemax-sitemin,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void GammaRateProcessVI::GlobalGetRateAlphaVI()  {

        if (! RateAlphaVI)  {
                 RateAlphaVI = new double[GetNsite()];
        }

        assert(GetMyid() == 0);
        int nprocs = GetNprocs();
        int i, width, smin[nprocs-1], smax[nprocs-1], workload[nprocs-1];
        MPI_Status stat;
        MESSAGE signal = RATEALPHAVI;
        
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        
        width = GetNsite()/(nprocs-1);
        for (i=0; i<nprocs-1; ++i) {
                 smin[i] = width*i;
                 smax[i] = width*(1+i);
                 if (i == (nprocs-2)) smax[i] = GetNsite();
        }
        for (i=1; i<nprocs; ++i)  {
                 MPI_Recv(RateAlphaVI+smin[i-1],smax[i-1]-smin[i-1],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        }
}

void GammaRateProcessVI::SlaveGetRateBetaVI()  {
 
        assert(GetMyid() > 0);
        
        int sitemin, sitemax, width = GetNsite()/(GetNprocs()-1);
        sitemin = (GetMyid()-1)*width;
        if (GetMyid() == (GetNprocs()-1))  {
                     sitemax = GetNsite();
        }
        else    {
                     sitemax = GetMyid()*width;
        }
        MPI_Send(RateBetaVI+sitemin,sitemax-sitemin,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void GammaRateProcessVI::GlobalGetRateBetaVI()  {

        if (! RateBetaVI)  {
                 RateBetaVI = new double[GetNsite()];
        }

        assert(GetMyid() == 0);
        int nprocs = GetNprocs();
        int i, width, smin[nprocs-1], smax[nprocs-1], workload[nprocs-1];
        MPI_Status stat;
        MESSAGE signal = RATEBETAVI;
        
        MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        
        width = GetNsite()/(nprocs-1);
        for (i=0; i<nprocs-1; ++i) {
                 smin[i] = width*i;
                 smax[i] = width*(1+i);
                 if (i == (nprocs-2)) smax[i] = GetNsite();
        }
        for (i=1; i<nprocs; ++i)  {
                 MPI_Recv(RateBetaVI+smin[i-1],smax[i-1]-smin[i-1],MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
        }
}


double GammaRateProcessVI::GlobalELBORate()  {

       assert(GetMyid() == 0);
       MPI_Status stat;
       MESSAGE signal = ELBORATEVI;
       MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
       
       int ratwidth = GetNsite()/(GetNprocs()-1);
       int ratmin[GetNprocs()-1];
       int ratmax[GetNprocs()-1];
       for (int i=0; i<GetNprocs()-1; ++i)  {
               ratmin[i] = ratwidth * i;
               ratmax[i] = ratwidth * (1+i);
               if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
       }
       ELBO_Rate = 0.0;
       double sum;
       for (int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(&sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
               ELBO_Rate += sum;
       }
       return ELBO_Rate;
}

void GammaRateProcessVI::SlaveELBORate()  {

      assert(GetMyid() > 0);
      
      int ratwidth = GetNsite()/(GetNprocs()-1);
      int ratmin[GetNprocs()-1];
      int ratmax[GetNprocs()-1];
      for (int i=0; i<GetNprocs()-1; ++i)  {
             ratmin[i] = ratwidth*i;
             ratmax[i] = ratwidth*(1+i);
             if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
      }
      
      double totalrate;
      for (int i=ratmin[GetMyid()-1]; i<ratmax[GetMyid()-1]; i++)  {

              totalrate += GetSiteRateSuffStatCount(i) * ((gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))));
              totalrate += (alpha * log(alpha) - rnd::GetRandom().logGamma(alpha) + (alpha - 1) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))) - alpha * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
              totalrate -= (GetRateAlphaVI(i) * log(GetRateBetaVI(i)) - rnd::GetRandom().logGamma(GetRateAlphaVI(i)));
              totalrate -= ((GetRateAlphaVI(i) - 1) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))) - GetRateBetaVI(i) * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
      }
      MPI_Send(&totalrate,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

/*void GammaRateProcessVI::SlaveELBORate()  {

      assert(GetMyid() > 0);
      
      int ratwidth = GetNsite()/(GetNprocs()-1);
      int ratmin[GetNprocs()-1];
      int ratmax[GetNprocs()-1];
      for (int i=0; i<GetNprocs()-1; ++i)  {
             ratmin[i] = ratwidth*i;
             ratmax[i] = ratwidth*(1+i);
             if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
      }
      double totratesuffstat = 0;
      for (int i=ratmin[GetMyid()-1]; i<ratmax[GetMyid()-1]; i++)     {
              totratesuffstat += GetSiteRateSuffStatCount(i);
      }
      double totalrate;
      for (int i=ratmin[GetMyid()-1]; i<ratmax[GetMyid()-1]; i++)  {
              double ratesuffstat = GetSiteRateSuffStatCount(i) / totratesuffstat;
              totalrate += ratesuffstat * ((gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))));
              totalrate += (alpha * log(alpha) - rnd::GetRandom().logGamma(alpha) + (alpha - 1) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))) - alpha * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
              totalrate -= (GetRateAlphaVI(i) * log(GetRateBetaVI(i)) - rnd::GetRandom().logGamma(GetRateAlphaVI(i)));
              totalrate -= ((GetRateAlphaVI(i) - 1) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))) - GetRateBetaVI(i) * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
      }
      MPI_Send(&totalrate,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}*/

void GammaRateProcessVI::SlaveGetMeanRateVI()  {

      assert(GetMyid() > 0);
 
      int ratwidth = GetNsite()/(GetNprocs()-1);
      int ratmin[GetNprocs()-1];
      int ratmax[GetNprocs()-1];
      for (int i=0; i<GetNprocs()-1; ++i)  {
              ratmin[i] = ratwidth*i;
              ratmax[i] = ratwidth*(i+1);
              if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
      } 

      double totalmeanrate;
      for (int j=ratmin[GetMyid()-1]; j<ratmax[GetMyid()-1]; j++)   {
              totalmeanrate += (GetRateAlphaVI(j)/GetRateBetaVI(j));
      }
      MPI_Send(&totalmeanrate,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

double GammaRateProcessVI::GlobalGetMeanRateVI()  {

      assert(GetMyid() == 0);
      MPI_Status stat;
      MESSAGE signal = MEANRATEVI;
      MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

      int ratwidth = GetNsite()/(GetNprocs()-1);
      int ratmin[GetNprocs()-1];
      int ratmax[GetNprocs()-1];
      for (int i=0; i<GetNprocs()-1; ++i)  {
              ratmin[i] = ratwidth*i;
              ratmax[i] = ratwidth*(1+i);
              if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
      }
      MeanRateVI = 0.0;
      double sum;
      for (int i=1; i<GetNprocs(); ++i)  {
              MPI_Recv(&sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
              MeanRateVI += sum;
      }
      return MeanRateVI;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void GammaRateProcessVI::GlobalMoveRate() {

        assert(GetMyid() == 0);
	MPI_Status stat;
	MESSAGE signal = UPDATE_MOVEPRRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        // Split Nsite among GetNprocs()-1 slaves
        int ratwidth = GetNsite()/(GetNprocs()-1);
        // int ratwidth = (GetSiteMax()-GetSiteMin())/(GetNprocs()-1);
        int ratmin[GetNprocs()-1];
        int ratmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i) {
                ratmin[i] = ratwidth*i;
                ratmax[i] = ratwidth*(1+i);
                if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
        }
        double* tmpratal = new double[GetNsite() + 1];
        for(int i=1; i<GetNprocs(); ++i) {
                MPI_Recv(tmpratal,(ratmax[i-1] - ratmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
                int k = 0;
                for(int j=ratmin[i-1]; j<ratmax[i-1]; ++j) {
                        rate[j] = tmpratal[k];
                        k++;
                }
        }
        delete[] tmpratal;
}

void GammaRateProcessVI::SlaveMoveRate() {

        assert(GetMyid() > 0);
        double* tmpratal = new double[GetNsite() + 1];
        
        // Split Nsite among GetNprocs()-1 slaves
        int ratwidth = GetNsite()/(GetNprocs()-1);
        // int ratwidth = (GetSiteMax()-GetSiteMin())/(GetNprocs()-1);
        int ratmin[GetNprocs()-1];
        int ratmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i) {
              ratmin[i] = ratwidth*i;
              ratmax[i] = ratwidth*(1+i);
              if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
        }
        int k = 0;
        for(int j=ratmin[GetMyid()-1]; j<ratmax[GetMyid()-1]; j++) {
              rate[j] = rnd::GetRandom().Gamma(GetRateAlphaVI(j), GetRateBetaVI(j));
              tmpratal[k] = rate[j];
              k++;  
        }
      
        MPI_Send(tmpratal,(ratmax[GetMyid()-1] - ratmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] tmpratal;
}

void GammaRateProcessVI::GlobalMoveRateParameter() {

        assert(GetMyid() == 0);
	MPI_Status stat;
	MESSAGE signal = UPDATE_MOVERATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        // Split Nsite among GetNprocs()-1 slaves
        int movwidth = GetNsite()/(GetNprocs()-1);
        // int movwidth = (GetSiteMax()-GetSiteMin())/(GetNprocs()-1);
        int movmin[GetNprocs()-1];
        int movmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i)  {
                movmin[i] = movwidth*i;
                movmax[i] = movwidth*(1+i);
                if (i == (GetNprocs()-2)) movmax[i] = GetNsite();
        }
        double* tmpmoval = new double[GetNsite() + 1];
        double* tmpmovbe = new double[GetNsite() + 1];
        for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpmoval,(movmax[i-1] - movmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int k = 0;
               for(int j=movmin[i-1]; j<movmax[i-1]; ++j) {
                       RateAlphaVI[j] = tmpmoval[k];
                       k++;
               }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpmovbe,(movmax[i-1] - movmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int l = 0;
               for(int j=movmin[i-1]; j<movmax[i-1]; ++j) {
                       RateBetaVI[j] = tmpmovbe[l];
                       l++;
               }
        }
        delete[] tmpmoval;
        delete[] tmpmovbe;
}

void GammaRateProcessVI::SlaveMoveRateParameter()  {
       
        assert(GetMyid() > 0);
       double* tmpmoval = new double[GetNsite() + 1];
       double* tmpmovbe = new double[GetNsite() + 1];

       // Split Nbranch among GetNprocs() - 1 slaves
       int movwidth = GetNsite()/(GetNprocs()-1);
       // int movwidth = (GetSiteMax()-GetSiteMin())/(GetNprocs()-1);
       int movmin[GetNprocs()-1];
       int movmax[GetNprocs()-1];
       for(int i=0; i<GetNprocs()-1; ++i) {
             movmin[i] = movwidth*i;
             movmax[i] = movwidth*(1+i);
             if (i == (GetNprocs()-2)) movmax[i] = GetNsite();
       }
       for(int j=movmin[GetMyid()-1]; j<movmax[GetMyid()-1]; j++) {
             LearningRateRate[j] = pow(taurate + j, -1 * kapparate); 
       }
       int k = 0;
       for(int j=movmin[GetMyid()-1]; j<movmax[GetMyid()-1]; j++) { 
              RateAlphaVI[j] = ((1 - LearningRateRate[j]) * RateAlphaVI[j]) + (LearningRateRate[j] * GetRateAlphahat(j)); 
              tmpmoval[k] = RateAlphaVI[j];
              k++;
       }
       int l = 0;
       for(int j=movmin[GetMyid()-1]; j<movmax[GetMyid()-1]; j++) {
              RateBetaVI[j] = ((1 - LearningRateRate[j]) * RateBetaVI[j]) + (LearningRateRate[j] * GetRateBetahat(j)); 
              tmpmovbe[l] = RateBetaVI[j];
              l++;
       } 

       MPI_Send(tmpmoval,(movmax[GetMyid()-1] - movmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       MPI_Barrier(MPI_COMM_WORLD);
       MPI_Send(tmpmovbe,(movmax[GetMyid()-1] - movmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       delete[] tmpmoval;
       delete[] tmpmovbe;
}

void GammaRateProcessVI::GlobalEstimateRateParameter()	{

	assert(GetMyid() == 0);
	MPI_Status stat;
	MESSAGE signal = UPDATE_ESTRATE;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        // Split Nsite among GetNprocs()-1 slaves
        int ratwidth = GetNsite()/(GetNprocs()-1);
        // int ratwidth = (GetSiteMax()-GetSiteMin())/(GetNprocs()-1);
        int ratmin[GetNprocs()-1];
        int ratmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i) {
                ratmin[i] = ratwidth*i;
                ratmax[i] = ratwidth*(1+i);
                if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
        }
        double* tmpratal = new double[GetNsite() + 1];
        double* tmpratbe = new double[GetNsite() + 1];
        for(int i=1; i<GetNprocs(); ++i) {
                MPI_Recv(tmpratal,(ratmax[i-1] - ratmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
                int k = 0;
                for(int j=ratmin[i-1]; j<ratmax[i-1]; ++j) {
                        RateAlphahat[j] = tmpratal[k];
                        k++;
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); ++i) {
                MPI_Recv(tmpratbe,(ratmax[i-1] - ratmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
                int l = 0;
                for(int j=ratmin[i-1]; j<ratmax[i-1]; ++j) {
                        RateBetahat[j] = tmpratbe[l];
                        l++;
                }
        }
        delete[] tmpratal;
        delete[] tmpratbe;  
}

void GammaRateProcessVI::SlaveEstimateRateParameter()	{

        assert(GetMyid() > 0);
        double* tmpratal = new double[GetNsite() + 1];
        double* tmpratbe = new double[GetNsite() + 1];
        
        // Split Nsite among GetNprocs()-1 slaves
        int ratwidth = GetNsite()/(GetNprocs()-1);
        // int ratwidth = (GetSiteMax()-GetSiteMin())/(GetNprocs()-1);
        int ratmin[GetNprocs()-1];
        int ratmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i) {
              ratmin[i] = ratwidth*i;
              ratmax[i] = ratwidth*(1+i);
              if (i == (GetNprocs()-2)) ratmax[i] = GetNsite();
        }
        int k = 0;
        for(int j=ratmin[GetMyid()-1]; j<ratmax[GetMyid()-1]; j++) {
              RateAlphahat[j] = alpha + GetSiteRateSuffStatCount(j);
              tmpratal[k] = RateAlphahat[j];
              k++;  
        }
        int l = 0;
        for(int j=ratmin[GetMyid()-1]; j<ratmax[GetMyid()-1]; j++) {
              RateBetahat[j] = beta + GetSiteRateSuffStatBeta(j);
              tmpratbe[l] = RateBetahat[j];
              l++;
        }   

        MPI_Send(tmpratal,(ratmax[GetMyid()-1] - ratmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Send(tmpratbe,(ratmax[GetMyid()-1] - ratmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        delete[] tmpratal;
        delete[] tmpratbe;
}

double GammaRateProcessVI::MoveRate()  {
        
        for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
        // for(int i=0; i<GetNsite(); i++) {
		rate[i] = rnd::GetRandom().Gamma(GetRateAlphaVI(i), GetRateBetaVI(i));
	}
        return 1.0;
}

double GammaRateProcessVI::MoveRateParameter()  {
       
       
       // for(int i=0; i<GetNsite(); i++) {
       for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{    
            LearningRateRate[i] = pow(taurate + i, -1 * kapparate);
            RateAlphaVI[i] = ((1 - LearningRateRate[i]) * RateAlphaVI[i]) + (LearningRateRate[i] * GetRateAlphahat(i));
            RateBetaVI[i] = ((1 - LearningRateRate[i]) * RateBetaVI[i]) + (LearningRateRate[i] * GetRateBetahat(i));
            // RateAlphaVI[i] = 1;
            // RateBetaVI[i] = 10; 
       }
   return 1.0;
}

void GammaRateProcessVI::EstimateRateParameter()  {
       
       // GlobalUpdateSiteRateSuffStat();
       // TotalParameterLength();
       for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
       // for(int i=0; i<GetNsite(); i++) {
            RateAlphahat[i] = alpha + GetSiteRateSuffStatCount(i);
            RateBetahat[i] = beta + GetSiteRateSuffStatBeta(i);
            // RateAlphahat[i] = 1;
            // RateBetahat[i] = 1;  
       }   
} 
