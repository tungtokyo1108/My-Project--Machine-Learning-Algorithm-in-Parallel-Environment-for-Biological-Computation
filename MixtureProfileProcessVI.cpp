
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#include "MixtureProfileProcessVI.h"
#include "Random.h"


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* MixtureProfileProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void MixtureProfileProcessVI::Create(int innsite, int indim)	{
	
       // PCAT = new double[GetNmodeMax()];
       if (! profile)	{
		ProfileProcess::Create(innsite,indim);
		allocprofile = new double[GetNmodeMax() * GetDim()];
		profile = new double*[GetNmodeMax()];
                dirweightVI = new double*[GetNmodeMax()];
                dirweighthat = new double*[GetNmodeMax()];
                dirweight_grad = new double*[GetNmodeMax()];
                LearningDirweight = new double[GetDim()];
		for (int i=0; i<GetNmodeMax(); i++)	{
			profile[i] = allocprofile + i*GetDim();
                        dirweightVI[i] = new double[GetDim()];
                        dirweighthat[i] = new double[GetDim()];
                        dirweight_grad[i] = new double[GetDim()]; 
			// profile[i] = new double[GetDim()];
		}
		alloc = new int[GetNsite()];
		occupancy = new int[GetNmodeMax()];
		dirweight = new double[GetDim()];
		logstatprior = new double[GetNmodeMax()];
		profilesuffstatlogprob = new double[GetNmodeMax()];
	}
}

void MixtureProfileProcessVI::Delete()	{
	if (profile)	{
		delete[] profilesuffstatlogprob;
		delete[] logstatprior;
		delete[] allocprofile;
		delete[] profile;
		profile = 0;
                for (int i=0; i<GetNmodeMax(); i++) {
                      delete[] dirweighthat[i];
                      delete[] dirweight_grad[i];
                }
                delete[] dirweighthat;
                delete[] dirweight_grad;
                // dirweighthat = 0;
                // dirweight_grad = 0;
		ProfileProcess::Delete();
	}
}

double MixtureProfileProcessVI::GetMeanDirWeight()	{
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
             for (int j=0; j<GetDim(); j++)  {
                  total += dirweightVI[i][j];
             }
        }
	return total;
}

double MixtureProfileProcessVI::GetStatEnt()	{
	double total = 0;
	UpdateOccupancyNumbers();
	for (int k=0; k<GetNcomponent(); k++)	{
		total += occupancy[k] * GetStatEnt(k);
	}
	return total / GetNsite();
}

double MixtureProfileProcessVI::GetStatEnt(int k)	{
	double total = 0;
	for (int i=0; i<GetDim(); i++)	{
		if (profile[k][i] <= 0)	{
			cerr << "error: 0 entry in profile\n";
			cerr << profile[k][i] << '\n';
			exit(1);
		}
		total -= profile[k][i] * log(profile[k][i]);
	}
	if (isnan(total))	{
		cerr << "entropy is nan\n";
		exit(1);
	}
	return  total;
}

double MixtureProfileProcessVI::GetCenterStatEnt()	{
	double totalweight = 0;
        for (int i=0; i<GetNmodeMax(); i++)	{
	         for (int k=0; k<GetDim(); k++)	{
		           totalweight += dirweightVI[i][k];
	         }
        }
	double total = 0;
        for (int i=0; i<GetNmodeMax(); i++)	{
 	         for (int k=0; k<GetDim(); k++)	{
		        double w = dirweightVI[i][k] / totalweight;
		        total -= w * log(w);
	         }
        }
	return total;
}

void MixtureProfileProcessVI::RenormalizeProfiles()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		double total = 0;
		for (int k=0; k<GetDim(); k++)	{
			total += profile[i][k];
		}
		for (int k=0; k<GetDim(); k++)	{
			profile[i][k] /= total;
		}
	}
}

void MixtureProfileProcessVI::SampleProfile()	{
	SampleHyper();
	SampleAlloc();
	SampleStat();
	// UpdateComponents();
}

void MixtureProfileProcessVI::SampleStat()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		SampleStat(profile[i]);
	}
}

void MixtureProfileProcessVI::SampleStat(int i)	{
	SampleStat(profile[i]);
}

double MixtureProfileProcessVI::ResampleEmptyProfiles()	{

	UpdateOccupancyNumbers();
	for (int i=0; i<GetNcomponent(); i++)	{
		if (! occupancy[i])	{
			SampleStat(i);
		}
	}
	return 1.0;
}
/*
void MixtureProfileProcess::SampleStat(double* prof)	{
	double total = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = rnd::GetRandom().sGamma(dirweight[k]);
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
	}
	total = 0;
	for (int k=0; k<GetDim(); k++)	{
		if (prof[k] < stateps)	{
			prof[k] = stateps;
		}
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
		if (isnan(prof[k]))	{
			cerr << "nan in sample stat\n";
			cerr << dirweight[k] << '\n';
			exit(1);
		}
	}
	// UpdateComponent(i);
}
*/
void MixtureProfileProcessVI::SampleStat(double* prof, double statmin)	{
	if (! statmin)	{
		statmin = stateps;
	}
	double total = 0;
	int infreached = 0;
	for (int k=0; k<GetDim(); k++)	{
		prof[k] = rnd::GetRandom().sGamma(dirweight[k]);
		if (prof[k] < statmin)	{
			prof[k] = statmin;
			infreached = 1;
		}
		total += prof[k];
	}
	for (int k=0; k<GetDim(); k++)	{
		prof[k] /= total;
	}
	if (infreached)	{
		statinfcount++;
	}
	totstatcount++;
}

double MixtureProfileProcessVI::LogProfilePrior()	{
	double total = 0;
	total += LogHyperPrior();
	total += LogAllocPrior();
	total += LogStatPrior();
	return total;
}

void MixtureProfileProcessVI::UpdateOccupancyNumbers()	{
	for (int i=0; i<GetNcomponent(); i++)	{
		occupancy[i] = 0;
	}
	for (int i=0; i<GetNsite(); i++)	{
		occupancy[alloc[i]]++;
	}
}

double MixtureProfileProcessVI::LogStatPrior()	{

	for (int i=0; i<GetNcomponent(); i++)	{
		LogStatPrior(i);
	}
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += logstatprior[i];
	}
	return total;
}

double MixtureProfileProcessVI::LogStatPrior(int cat)	{
	double total = 0;
	double totalweight = 0;
	for (int k=0; k<GetDim(); k++)	{
		total += (dirweight[k] - 1) * log(profile[cat][k]) - rnd::GetRandom().logGamma(dirweight[k]);
		totalweight += dirweight[k];
	}
	total += rnd::GetRandom().logGamma(totalweight);
	logstatprior[cat] = total;
	return total;
}

double MixtureProfileProcessVI::ProfileSuffStatLogProb()	{
	// simply, sum over all components
	for (int i=0; i<GetNcomponent(); i++)	{
		ProfileSuffStatLogProb(i);
	}
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		total += profilesuffstatlogprob[i];
	}
	return total;
}

double MixtureProfileProcessVI::MoveDirWeights(double tuning, int nrep)	{
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		for (int k=0; k<GetDim(); k++)	{
			double deltalogprob = - LogHyperPrior() - LogStatPrior();
			double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
			double e = exp(m);
			dirweight[k] *= e;
			deltalogprob += LogHyperPrior() + LogStatPrior();
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
	ResampleEmptyProfiles();
	return naccepted / nrep / GetDim();
}

void MixtureProfileProcessVI::SwapComponents(int cat1, int cat2)	{

	int tmp = occupancy[cat1];
	occupancy[cat1] = occupancy[cat2];
	occupancy[cat2] = tmp;

	for (int k=0; k<GetDim(); k++)	{
		double tmp = profile[cat1][k];
		profile[cat1][k] = profile[cat2][k];
		profile[cat2][k] = tmp;
	}
	/*
	double* temp = profile[cat1];
	profile[cat1] = profile[cat2];
	profile[cat2] = temp;
	*/

	for (int i=0; i<GetNsite(); i++)	{
		if (alloc[i] == cat1)	{
			alloc[i] = cat2;
		}
		else if (alloc[i] == cat2)	{
			alloc[i] = cat1;
		}
	}
}

