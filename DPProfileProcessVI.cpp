
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#include "DPProfileProcessVI.h"
#include "Random.h"
#include <gsl/gsl_sf_psi.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* DPProfileProcessVI
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

const double meandir[] = {0.499737,0.171262,0.183399,0.225593,0.197453,0.211819,0.173191,0.175454,0.3181,0.240008,0.187577,0.324778,0.205587,0.395097,0.162356,0.519427,0.526213,0.349177,0.0511527,0.130222};

void DPProfileProcessVI::SampleHyper()	{
	kappa = 1;

        /*kappa_alpha = new double[GetNmodeMax()];
        kappa_beta = new double[GetNmodeMax()];
        for (int k=0; k<GetNmodeMax(); k++)	{
		kappa_alpha[k] = 1;
                kappa_beta[k] = 1;
	}*/
	// kappa = GetNsite() / 5 ;
	for (int i=0; i<GetDim(); i++)	{
		dirweight[i] = 1;
		// dirweight[i] = meandir[i];
                LearningDirweight[i] = 0;
	}
        for (int i=0; i<GetNmodeMax(); i++)	{
             for (int j=0; j<GetDim(); j++)  {
                  dirweightVI[i][j] = rnd::GetRandom().sExpo();
             }
        }
        /*for (int i=0; i<GetNmodeMax(); i++)	{
             for (int j=0; j<GetDim(); j++)  {
                  dirweight_grad[i][j] = 1e-8;
             }
        }*/
        /*double tot = 0;
            for (int i=0; i<GetNcomponent(); i++)  {
                  for(int j=0; j<GetDim(); j++)        {
                       tot += dirweightVI[i][j];
                  }
                  double total = 0; 
                  for (int j=0; j<GetDim(); j++)       {
                       total += gsl_sf_psi(dirweightVI[i][j]) - gsl_sf_psi(tot); 
                  }
                  PCAT[i] = exp(total);
            }*/
}
	
// 1 component
/*
void DPProfileProcess::SampleAlloc()	{
	CreateComponent(0);
	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{
		AddSite(i,0);
	}
}
*/
	
// Nsite compoents
/*
void DPProfileProcess::SampleAlloc()	{
	Ncomponent = GetNsite();
	for (int i=0; i<GetNsite(); i++)	{
		CreateComponent(i);
		AddSite(i,i);
	}
}
*/

/*double DPProfileProcessVI::LogProxy(int site, int cat)	{
	return 0;
}*/

void DPProfileProcessVI::SampleAlloc()	{

	CreateComponent(0);
	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{
		double* p = new double[Ncomponent+1];
		double total = 0;
		double max = 0;
		for (int k=0; k<Ncomponent; k++)	{
			double tmp = log(occupancy[k]) * LogProxy(i,k);
			if ((!k) || (max < tmp))	{
				max = tmp;
			}
			p[k] = tmp;
		}
		p[Ncomponent] = log(kappa) + LogProxy(i,Ncomponent);
		if (max < p[Ncomponent])	{
			max = p[Ncomponent];
		}
		for (int k=0; k<=Ncomponent; k++)	{
			double tmp = exp(p[k] - max);
			total += tmp;
			p[k] = total;
		}
		double q = total * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<=Ncomponent) && (q > p[k])) k++;
		if (k == Ncomponent+1)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		if (k==Ncomponent)	{
			CreateComponent(k);
			Ncomponent++;
		}
		AddSite(i,k);
		delete[] p;
	}
}

double DPProfileProcessVI::LogHyperPrior()	{
	double total = 0;
	if (kappaprior == 0)	{
		total = -kappa / 10.0;
	}
	else 	{
		total = -log(kappa);
		if ((kappa < 1e-4) || (kappa > 1e4))	{
			// total = InfProb;
			total -= 1.0 / 0;
		}
	}
	double sum = 0;
	for (int k=0; k<GetDim(); k++)	{
		total -= dirweight[k];
		sum += dirweight[k];
	}
	if (sum < GetMinTotWeight())	{
		// total += InfProb;
		total -= 1.0 / 0;
	}
	return total;
}

double DPProfileProcessVI::LogAllocPrior()	{
	double total = GetNOccupiedComponent() * log(kappa);
	// double total = GetNcomponent() * log(kappa);
	for (int i=0; i<GetNsite(); i++)	{
		total -= log(kappa + i);
	}
	return total;
}

double DPProfileProcessVI::MoveHyper(double tuning, int nrep)	{
	double total = 0;
	if (!burnin)	{
		total += MoveKappa(tuning,nrep);
	}
	total += MoveDirWeights(tuning,nrep);
	return total;
}

double DPProfileProcessVI::MoveKappa(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogAllocPrior();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		kappa *= e;
		deltalogprob += LogHyperPrior() + LogAllocPrior();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		/*
		if (kappa < 1)	{
			accepted = 0;
		}
		*/
		if (accepted)	{
			naccepted++;
		}
		else	{
			kappa /= e;
		}
	}
	return naccepted / nrep;
}



