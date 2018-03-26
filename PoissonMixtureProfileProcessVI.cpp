
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#include "PoissonMixtureProfileProcessVI.h"
#include "Random.h"
#include <gsl/gsl_sf_psi.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PoissonMixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void PoissonMixtureProfileProcessVI::Create(int innsite, int indim)	{

       /*dirweightVI = new double*[GetNmodeMax()];
       dirweighthat = new double*[GetNmodeMax()];
       for (int i=0; i<GetNmodeMax(); i++)  {
             dirweightVI[i] = new double[GetDim()];
             dirweighthat[i] = new double[GetDim()];
       }

       for (int i=0; i<GetNmodeMax(); i++)  {
             for (int j=0; j<GetDim(); j++)     {
                    dirweightVI[i][j] = 1;  
             } 
       }*/
       if (! profilesuffstatcount)	{
		PoissonProfileProcess::Create(innsite,indim);
		MixtureProfileProcessVI::Create(innsite,indim);
		profilesuffstatcount  = new int*[GetNmodeMax()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profilesuffstatcount[i] = new int[GetDim()];
		}
		// SampleProfile();
	}
}

void PoissonMixtureProfileProcessVI::Delete() {
	if (profilesuffstatcount)	{
		for (int i=0; i<GetNmodeMax(); i++)	{
			delete[] profilesuffstatcount[i];
		}
		delete[] profilesuffstatcount;
		profilesuffstatcount = 0;
		PoissonProfileProcess::Delete();
		MixtureProfileProcessVI::Delete();
	}
}

void PoissonMixtureProfileProcessVI::UpdateModeProfileSuffStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[i][k] = 0;
		}
	}
	for (int i=0; i<GetNsite(); i++)	{
		const int* count = GetSiteProfileSuffStatCount(i);
		int cat = alloc[i];
		for (int k=0; k<GetDim(); k++)	{
			profilesuffstatcount[cat][k] += count[k];
		}
	}
}

double PoissonMixtureProfileProcessVI::DiffLogSampling(int cat, int site)	{

	const int* nsub = GetSiteProfileSuffStatCount(site);
	int* catnsub = profilesuffstatcount[cat];
	int totalsub = 0;
	double priorweight = 0;
	int grandtotal = 0;
	for (int k=0; k<GetDim(); k++)	{
		totalsub += nsub[k];
		priorweight += dirweight[k];
		grandtotal += catnsub[k];
	}
	
	double total = 0;
	for (int j=0; j< totalsub; j++)	{
		total -= log(priorweight + grandtotal + j);
	}
	for (int k=0; k<GetDim(); k++)	{
		for (int j=0; j< nsub[k]; j++)	{
			total += log(dirweight[k] + catnsub[k] + j);
		}
	}
	return total;
}

double PoissonMixtureProfileProcessVI::LogStatProb(int site, int cat)	{
	
       const int* nsub = GetSiteProfileSuffStatCount(site);
       double total = 0;
       double tot = 0; 
                for (int k=0; k<GetDim(); k++)  {
                      tot += dirweightVI[cat][k];  
                }
		for (int k=0; k<GetDim(); k++)	{
		      total += nsub[k] * (gsl_sf_psi(dirweightVI[cat][k]) - gsl_sf_psi(tot));
                      // total -= (dirweightVI[cat][k] / tot);
		}
	return total;
}

/*double PoissonMixtureProfileProcessVI::LogStatProb(int site, int cat)	{
	
       const int* nsub = GetSiteProfileSuffStatCount(site);
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += nsub[k] * log(profile[cat][k]);
	}
	return total;
}*/

void PoissonMixtureProfileProcessVI::RemoveSite(int site, int cat)	{
	occupancy[cat] --;
	if (activesuffstat)	{
		const int* nsub = GetSiteProfileSuffStatCount(site);
		int* catnsub = profilesuffstatcount[cat];
		for (int k=0; k<GetDim(); k++)	{
			catnsub[k] -= nsub[k];
		}
	}
}

void PoissonMixtureProfileProcessVI::AddSite(int site, int cat)	{
	alloc[site] = cat;
	occupancy[cat] ++;
	UpdateZip(site);
	if (activesuffstat)	{
		const int* nsub = GetSiteProfileSuffStatCount(site);
		int* catnsub = profilesuffstatcount[cat];
		for (int k=0; k<GetDim(); k++)	{
			catnsub[k] += nsub[k];
		}
	}
}

void PoissonMixtureProfileProcessVI::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcessVI::SwapComponents(cat1,cat2);
	for (int k=0; k<GetDim(); k++)	{
		int tmp = profilesuffstatcount[cat1][k];
		profilesuffstatcount[cat1][k] = profilesuffstatcount[cat2][k];
		profilesuffstatcount[cat2][k] = tmp;
	}
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                                     Variational distribution 
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double PoissonMixtureProfileProcessVI::ProfileSuffStatLogProb(int cat)	{
	double total = 0;
	double priorweight = 0;
	double postweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[k] + profilesuffstatcount[cat][k]) - rnd::GetRandom().logGamma(dirweight[k]);
		priorweight += dirweight[k];
		postweight += dirweight[k] + profilesuffstatcount[cat][k];
	}
	total += rnd::GetRandom().logGamma(priorweight) - rnd::GetRandom().logGamma(postweight);
	return total;
}

double PoissonMixtureProfileProcessVI::TotDirWeightVI(int cat, int dim)  {

       double totVI = 0;
       for (int a = 0; a<cat; a++)  {
            for (int b = 0; b<dim; b++)  {
                  totVI += gsl_sf_psi(dirweightVI[a][b]);
            }
       }
       return totVI;
}

/*double PoissonMixtureProfileProcessVI::ELBODirWeight()  {
        
       GetTotalSitePCAT();
       ELBO_DirWeight = 0;
       double priorweight = 0;
       double priorweightVI = 0;
       double tot = 0;
       double totvariational = 0;
       for (int k=0; k<GetNcomponent(); k++)	{
                for (int i=0; i<GetDim(); i++)    {
                     tot += dirweightVI[k][i];
                }
                for (int i=0; i<GetDim(); i++)    {
                     ELBO_DirWeight += ((GetPCAT(k) * profilesuffstatcount[k][i]) * (gsl_sf_psi(dirweightVI[k][i]) - gsl_sf_psi(tot)));
                     ELBO_DirWeight += ((dirweight[i] -1) * (gsl_sf_psi(dirweightVI[k][i]) - gsl_sf_psi(tot)) - rnd::GetRandom().logGamma(dirweight[i]));
                     priorweight += dirweight[i];
                }
                ELBO_DirWeight += rnd::GetRandom().logGamma(priorweight);
       }
       for (int k=0; k<GetNcomponent(); k++)     {
                for (int i=0; i<GetDim(); i++)   {
                     totvariational += dirweightVI[k][i];
                }
                for (int i=0; i<GetDim(); i++)       {
                     ELBO_DirWeight -= (dirweightVI[k][i] - 1) * (gsl_sf_psi(dirweightVI[k][i]) - gsl_sf_psi(totvariational)) - rnd::GetRandom().logGamma(dirweightVI[k][i]);
                     priorweightVI += dirweightVI[k][i];
                }
                ELBO_DirWeight -= rnd::GetRandom().logGamma(priorweightVI);
       }
       return ELBO_DirWeight;
}*/

double PoissonMixtureProfileProcessVI::ELBODirWeight()   {

       ELBO_DirWeight = 0;
       for (int k=0; k<GetNcomponent(); k++)     {
               if (occupancy[k])    {
                      ELBO_DirWeight += ELBODirWeight(k);
               }
       }
       ELBO_DirWeight += PCATDirWeight();
       return ELBO_DirWeight;
}

double PoissonMixtureProfileProcessVI::ELBODirWeight(int cat)   {

       double tot_ELBO_DirWeight = 0;
       double priorweight = 0;
       double priorweightVI = 0;
       double tot = 0;
       double totvariational = 0;
       for (int k=0; k<GetDim(); k++)    {
               tot += dirweightVI[cat][k];
       }
       for (int k=0; k<GetDim(); k++)    {
               // tot_ELBO_DirWeight += ( (GetPCAT(cat) * profilesuffstatcount[cat][k]) * (gsl_sf_psi(dirweightVI[cat][k]) - gsl_sf_psi(tot)) );
               tot_ELBO_DirWeight += ( (dirweight[k] - 1) * (gsl_sf_psi(dirweightVI[cat][k]) - gsl_sf_psi(tot)) - rnd::GetRandom().logGamma(dirweight[k]));
               priorweight += dirweight[k];
       }
       tot_ELBO_DirWeight += rnd::GetRandom().logGamma(priorweight);
       for (int k=0; k<GetDim(); k++)    {
               tot_ELBO_DirWeight -= ( (dirweightVI[cat][k] - 1) * (gsl_sf_psi(dirweightVI[cat][k]) - gsl_sf_psi(tot)) - rnd::GetRandom().logGamma(dirweightVI[cat][k]) );
               priorweightVI += dirweightVI[cat][k];
       }
       tot_ELBO_DirWeight -= rnd::GetRandom().logGamma(priorweightVI);
       return tot_ELBO_DirWeight;
}

double PoissonMixtureProfileProcessVI::PCATDirWeight()   {

       double tot = 0;
       for (int k=0; k<GetNcomponent(); k++)     {
               if (occupancy[k])    {
                      tot += GetPCAT(k) * PCATDirWeight(k);
               }
       }
       return tot;
}

double PoissonMixtureProfileProcessVI::PCATDirWeight(int cat)    {

       double total = 0;
       double tot = 0;
       for (int k=0; k<GetDim(); k++)    {
               tot += dirweightVI[cat][k];
       }
       for (int k=0; k<GetDim(); k++)    {
               total += profilesuffstatcount[cat][k] * ( gsl_sf_psi(dirweightVI[cat][k]) - gsl_sf_psi(tot) );
       }
       return total;
}

/*double PoissonMixtureProfileProcessVI::ELBODirWeight()   {

       GetTotalSitePCAT();
       ELBO_DirWeight = 0;
       double priorweight = 0;
       double priorweightVI = 0;
       double tot = 0;
       double totvariational = 0;
       double totprofilesuffstat = 0;
       for (int k=0; k<GetNcomponent(); k++)   {
            for (int i=0; i<GetDim(); i++)        {
                 totprofilesuffstat += profilesuffstatcount[k][i] * GetPCAT(k);
            }  
       }
       for (int k=0; k<GetNcomponent(); k++)   {
            for (int i=0; i<GetDim(); i++)         {
                 tot += dirweightVI[k][i];
            }
            for (int i=0; i<GetDim(); i++)         {
                 double profilesuffstat = (profilesuffstatcount[k][i] * GetPCAT(k)) / totprofilesuffstat;
                 ELBO_DirWeight += ((profilesuffstat) * (gsl_sf_psi(dirweightVI[k][i]) - gsl_sf_psi(tot)));
                 ELBO_DirWeight += (dirweight[i] - 1) * (gsl_sf_psi(dirweightVI[k][i]) - gsl_sf_psi(tot)) - rnd::GetRandom().logGamma(dirweight[i]);
                 priorweight += dirweight[i];
            }
       }
       ELBO_DirWeight += rnd::GetRandom().logGamma(priorweight);
       for (int k=0; k<GetNcomponent(); k++)   {
            for (int i=0; i<GetDim(); i++)        {
                 totvariational += dirweightVI[k][i];
            }
            for (int i=0; i<GetDim(); i++)        {
                 ELBO_DirWeight -= (dirweightVI[k][i] - 1) * (gsl_sf_psi(dirweightVI[k][i]) - gsl_sf_psi(totvariational)) - rnd::GetRandom().logGamma(dirweightVI[k][i]);
                 priorweightVI += dirweightVI[k][i];
            }
       }     
       ELBO_DirWeight -= rnd::GetRandom().logGamma(priorweightVI);
       return ELBO_DirWeight;
}*/


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                            Stochastic Optimization for Variational Parameter
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double PoissonMixtureProfileProcessVI::GetProfileParameter(int site, int cat)  {

       const int* nsub = GetSiteProfileSuffStatCount(site);
       double tot = 0;
       double tot2 = 0; 
                for (int i=0; i<GetDim(); i++)  {
                      tot2 += dirweightVI[cat][i];  
                }
		for (int i=0; i<GetDim(); i++)	{
		      tot += nsub[i] * (gsl_sf_psi(dirweightVI[cat][i]) - gsl_sf_psi(tot2));
		}
       return tot;
}

double PoissonMixtureProfileProcessVI::Getprofilesuffstatcount(int cat) {
        UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();
        double tot = 0;
          for (int k=0; k<GetDim(); k++)  {
             tot += profilesuffstatcount[cat][k];
          }
    return tot;
}

double PoissonMixtureProfileProcessVI::EstimateDirWeight() {
        
        UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();
        GetTotalSitePCAT();
        double totweight = 0;
        
        // for (int k=0; k<GetNmodeMax(); k++)	{
        for (int k=0; k<GetNcomponent(); k++)	{
                for (int i=0; i<GetDim(); i++)    {
                       dirweighthat[k][i] = dirweight[i] + (profilesuffstatcount[k][i] * GetPCAT(k));
                       dirweight_grad[k][i] = dirweighthat[k][i] * dirweighthat[k][i];
                       // dirweighthat[k][i] = 10;
                }
	}
        // for (int k=0; k<GetNmodeMax(); k++)	{
        for (int k=0; k<GetNcomponent(); k++)	{
                for (int i=0; i<GetDim(); i++)    {
                       if (k==i) {
                             LearningDirweight[i] += dirweight_grad[k][i];   
                       }
                }
        }
   return 1.0;
}

double PoissonMixtureProfileProcessVI::MoveDirWeightVI()  {

       EstimateDirWeight();
       /*double taudirweight = rnd::GetRandom().sExpo();
       double kappadirweight = rnd::GetRandom().sExpo();
       double* LearningDirWeight = new double[GetDim()];*/
     // for (int k=0; k<GetNmodeMax(); k++)  {
     for (int k=0; k<GetNcomponent(); k++)	{
       for (int i=0; i<GetDim(); i++)        {
            /*LearningDirWeight[i] = pow(taudirweight + i, -1 * kappadirweight);
            dirweightVI[k][i] = (1 - LearningDirWeight[k]) * dirweightVI[k][i];
            dirweightVI[k][i] += LearningDirWeight[k] * dirweighthat[k][i];*/
            // dirweight_grad[i] += dirweighthat[k][i] * dirweighthat[k][i];
            dirweightVI[k][i] += dirweighthat[k][i] / (sqrt(LearningDirweight[i] + 1e-8));
       }
     }
     ELBODirWeight();
   return 1.0;
}

double PoissonMixtureProfileProcessVI::MoveProfile()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		MoveProfile(i);
	}
	return 1;
}

double PoissonMixtureProfileProcessVI::MoveProfile(int cat)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		profile[cat][k] = rnd::GetRandom().sGamma(dirweightVI[cat][k]);
		if (profile[cat][k] < stateps)	{
			profile[cat][k] = stateps;
		}
		total += profile[cat][k];
	}
	for (int k=0; k<GetDim(); k++)	{
		profile[cat][k] /= total;
	}
	return 1;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/* 


double PoissonMixtureProfileProcess::LogStatIntPrior()	{
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		if (occupancy[k])	{
			total += LogStatIntPrior(k);
		}
	}
	return total;
}

double PoissonMixtureProfileProcess::LogStatIntPrior(int cat)	{

	double total = 0;
	int tot = 0;
	double totweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += rnd::GetRandom().logGamma(dirweight[k] + profilesuffstatcount[cat][k]);
		total -= rnd::GetRandom().logGamma(dirweight[k]);
		totweight += dirweight[k];
		tot += profilesuffstatcount[cat][k];
	}
	total -= rnd::GetRandom().logGamma(totweight + tot);
	total += rnd::GetRandom().logGamma(totweight);
	return total;
}
	
double PoissonMixtureProfileProcess::MoveDirWeights(double tuning, int nrep)	{

	UpdateOccupancyNumbers();
	UpdateModeProfileSuffStat();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatIntPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[k] *= e;
			deltalogprob += LogHyperPrior() + LogStatIntPrior();
			deltalogprob += m;
			int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
			if (accepted)	{
				naccepted++;
			}
			else	{
				dirweight[k] /= e;
			}
		}
	}
	SampleStat();
	return naccepted / nrep / GetDim();
}

*/
