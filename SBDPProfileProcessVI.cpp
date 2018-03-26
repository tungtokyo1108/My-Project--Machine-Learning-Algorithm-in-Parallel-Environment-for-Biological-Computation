
/********************


**********************/


#include "SBDPProfileProcessVI.h"
#include "Random.h"
#include <gsl/gsl_sf_psi.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* SBDPProfileProcessVI
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


void SBDPProfileProcessVI::Create(int innsite, int indim)	{
        
        kappa_alpha = new double[GetNmodeMax()];
        kappa_beta = new double[GetNmodeMax()];
        
        for (int k=0; k<GetNmodeMax(); k++)	{
		kappa_alpha[k] = 1;
                kappa_beta[k] = GetNsite() / 5;
	}

	if (! V)	{
		DPProfileProcessVI::Create(innsite,indim);
		V = new double[GetNmodeMax()];
                V_P = new double[GetNmodeMax()];
		weight = new double[GetNmodeMax()];
                weightVI = new double[GetNmodeMax()];
	}
}

void SBDPProfileProcessVI::Delete()	{

        delete[] kappa_alpha;
        delete[] kappa_beta;
        // kappa_alpha = 0;
        // kappa_beta = 0; 
	if (V)	{
		delete[] V;
		delete[] weight;
                delete[] V_P;
                delete[] weightVI;
                // delete[] weight_P;
		DPProfileProcessVI::Delete();
	}
}

void SBDPProfileProcessVI::SampleAlloc()	{

	for (int k=0; k<GetNmodeMax(); k++)	{
		CreateComponent(k);
	}
	Ncomponent = GetNmodeMax();

	SampleWeights();
	for (int i=0; i<GetNsite(); i++)	{
		double U = rnd::GetRandom().Uniform();
		double total = weight[0];
		int k = 0;
		while ((k<GetNmodeMax()) && (total < U))	{
			k++;
			total += weight[k];
		}
		if (k == GetNmodeMax())	{
			cerr << "error in SBDPProfileProcess::SampleAlloc: overflow\n";
			exit(1);
		}
		AddSite(i,k);
	}
}


void SBDPProfileProcessVI::DrawProfileFromPrior()	{

	if (! GetMyid())	{
		cerr << "error: in master DrawProfileFromPrior\n";
		exit(1);
	}

	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
		RemoveSite(i,alloc[i]);
		int choose = rnd::GetRandom().FiniteDiscrete(GetNcomponent(),weight);
		AddSite(i,choose);
	}
}

void SBDPProfileProcessVI::IncrementalSampleAlloc()	{

	kappa = 0.1;

	for (int i=0; i<GetNsite(); i++)	{
		RemoveSite(i,alloc[i]);
	}

	AddSite(0,0);
	Ncomponent = 1;
	
	for (int i=0; i<GetNsite(); i++)	{

		int K = Ncomponent + 1;
		if (K > GetNmodeMax())	{
			K--;
		}
		double* p = new double[K];
		double total = 0;
		double max = 0;
		for (int k=0; k<K; k++)	{
			double w = occupancy[k];
			if (! w)	{
				w = kappa;
			}
			double tmp = log(w) * LogProxy(i,k);
			if ((!k) || (max < tmp))	{
				max = tmp;
			}
			p[k] = tmp;
		}
		for (int k=0; k<K; k++)	{
			double tmp = exp(p[k] - max);
			total += tmp;
			p[k] = total;
		}
		double q = total * rnd::GetRandom().Uniform();
		int k = 0;
		while ((k<K) && (q > p[k])) k++;
		if (k == K)	{
			cerr << "error in draw dp mode: overflow\n";
			exit(1);
		}
		if (k==Ncomponent)	{
			if (Ncomponent <= GetNmodeMax())	{
				Ncomponent++;
			}
		}
		AddSite(i,k);
		delete[] p;
	}

	Ncomponent = GetNmodeMax();
	ResampleWeights();
	cerr << "init incremental ok\n";
}

double SBDPProfileProcessVI::LogStatPrior()	{

	UpdateOccupancyNumbers();
	double total = 0;
	for (int i=0; i<GetNcomponent(); i++)	{
		if (occupancy[i])	{
			total += DPProfileProcessVI::LogStatPrior(i);
		}
	}
	return total;
}

void SBDPProfileProcessVI::SwapComponents(int cat1, int cat2)	{

	MixtureProfileProcessVI::SwapComponents(cat1,cat2);
	double tempv = V[cat1];
	V[cat1] = V[cat2];
	V[cat2] = tempv;
	double tempw = weight[cat1];
	weight[cat1] = weight[cat2];
	weight[cat2] = tempw;
}

void SBDPProfileProcessVI::SampleWeights()	{
        
        
	double cumulProduct = 1.0;
	double totweight = 0;
	double v, x, y;
	for (int k=0; k<GetNcomponent(); k++)	{
		x = rnd::GetRandom().sGamma(1.0);
		y = rnd::GetRandom().sGamma(kappa);
		v = x / (x+y);
		V[k] = v;
		if (k == GetNcomponent() - 1)	{
			V[k] = 1;
			v = 1;
		}
		weight[k] = v * cumulProduct;
		cumulProduct *= (1 - v);	
		totweight += weight[k];
	}
}


void SBDPProfileProcessVI::ResampleWeights()	{

	UpdateOccupancyNumbers();
	// ???
	// int remainingOcc = GetNsite();
	double cumulProduct = 0;
	// double totweight = 0;
	double v;
	for (int k=0; k<GetNcomponent(); k++)	{
   
            v = (gsl_sf_psi(kappa_alpha[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k]));
            weightVI[k] = v + cumulProduct;
            cumulProduct += (gsl_sf_psi(kappa_beta[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k]));	
	}
}

double SBDPProfileProcessVI::ELBOWeight()  {
 
       // GetTotalSitePCAT();
       ELBO_Weight = 0;
       
       double v;
       for (int k=0; k<GetNcomponent(); k++)  {
            v = GetPCAT(k) * (gsl_sf_psi(kappa_alpha[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k]));
            double cumulProduct = 0;
            for (int j=k+1; j<GetNcomponent(); j++)  {
                    cumulProduct += GetPCAT(j);
            }
            cumulProduct *= (gsl_sf_psi(kappa_beta[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k]));
            ELBO_Weight = v + cumulProduct; 
            ELBO_Weight -= GetPCAT(k) * log(GetPCAT(k));
       }
       return ELBO_Weight;
}

double SBDPProfileProcessVI::ELBOkappa()  {

       ELBO_kappa = 0;
       for (int k=0; k<GetNcomponent(); k++)  {
            ELBO_kappa += (rnd::GetRandom().logGamma(1 + kappa) - rnd::GetRandom().logGamma(1.0) - rnd::GetRandom().logGamma(kappa) + (1.0 - 1) * (gsl_sf_psi(kappa_alpha[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k])));
            ELBO_kappa += ((kappa - 1) * (gsl_sf_psi(kappa_beta[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k])));
            ELBO_kappa -= (rnd::GetRandom().logGamma(kappa_alpha[k] + kappa_beta[k]) - rnd::GetRandom().logGamma(kappa_alpha[k]) - rnd::GetRandom().logGamma(kappa_beta[k]));
            ELBO_kappa -= ((kappa_alpha[k] - 1) * (gsl_sf_psi(kappa_alpha[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k])) + (kappa_beta[k] - 1) * (gsl_sf_psi(kappa_beta[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k])));
       }
       return ELBO_kappa;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                                    Stochastic Optimization for Variational Parameter 
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double SBDPProfileProcessVI::GetSBDPParameter(int k) {
       
       double total = 0;
                for (int i=0; i<k-1; i++) {
                      total += gsl_sf_psi(kappa_beta[i]) - gsl_sf_psi(kappa_alpha[i] + kappa_beta[i]);
                }
                total += gsl_sf_psi(kappa_alpha[k]) - gsl_sf_psi(kappa_alpha[k] + kappa_beta[k]);
       return total; 
}

double SBDPProfileProcessVI::MoveKappa(double tuning, int nrep) {
   
   GetTotalSitePCAT();
   for(int rep=0; rep<nrep; rep++) {
       for (int k=0; k<GetNcomponent(); k++)	{
             kappa_alpha[k] = 1 + GetPCAT(k);
             // kappa_alpha[k] = 1;
             double total = 0;
                for(int j=k+1; j<GetNcomponent(); j++) {
                    total += GetPCAT(j);
                    // total += 0;  
                } 
             kappa_beta[k] = kappa + total;  
             // kappa_beta[k] = kappa;           
       }
    }
    return 1.0;   
}


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*

double SBDPProfileProcess::LogIntegratedAllocProb()	{
	int remainingOcc = GetNsite();
	double total = 0;
	for (int k=0; k<GetNcomponent(); k++)	{
		if (remainingOcc)	{
			remainingOcc -= occupancy[k];
			total += log(kappa) + rnd::GetRandom().logGamma(1 + occupancy[k]) + rnd::GetRandom().logGamma(kappa + remainingOcc) - rnd::GetRandom().logGamma(1 + kappa + occupancy[k] + remainingOcc);
		}
	}
	if (remainingOcc)	{
		cerr << "error in allocation count\n";
		exit(1);
	}
	return total;
}

double SBDPProfileProcess::MoveKappa(double tuning, int nrep)	{
	UpdateOccupancyNumbers();
	double naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogIntegratedAllocProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		kappa *= e;
		deltalogprob += LogHyperPrior() + LogIntegratedAllocProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			naccepted++;
		}
		else	{
			kappa /= e;
		}
	}
	// not useful
	// ResampleWeights();
	return naccepted / nrep;
}


double SBDPProfileProcess::MoveOccupiedCompAlloc(int k0)	{

	int nrep = (int) (k0 * kappa);
	UpdateOccupancyNumbers();
	ResampleWeights();
	double total = 0.0;
	int Nocc = GetNOccupiedComponent();
	if (Nocc != 1)	{
		for (int i=0; i<nrep; i++)	{
			int* occupiedComponentIndices = new int[Nocc];
			int j=0;
			for (int k=0; k<GetNcomponent(); k++)	{
				if (occupancy[k] != 0)	{
					occupiedComponentIndices[j] = k;
					j++;
				}
			}
			if (j != Nocc)	{
				cerr << "error in MoveOccupiedCompAlloc.\n";
				exit(1);
			}
			int* indices = new int[2];
			rnd::GetRandom().DrawFromUrn(indices,2,Nocc);
			int cat1 = occupiedComponentIndices[indices[0]];
			int cat2 = occupiedComponentIndices[indices[1]];
			double logMetropolis = (occupancy[cat2] - occupancy[cat1]) * log(weight[cat1] / weight[cat2]);
			int accepted = (log(rnd::GetRandom().Uniform()) < logMetropolis);
			if (accepted)	{
				total += 1.0;
				// SwapComponents(cat1, cat2);
				MixtureProfileProcess::SwapComponents(cat1,cat2);
			
			}
			delete[] occupiedComponentIndices;
			delete[] indices; 
		}
		return total /= nrep;
	}
	return 0;
}


double SBDPProfileProcess::MoveAdjacentCompAlloc(int k0)	{

	int nrep = (int) (k0 * kappa);
	ResampleWeights();
	
	double total = 0;

	for (int i=0; i<nrep; i++)	{
		//int cat1 = (int)(rnd::GetRandom().Uniform() * (GetNcomponent()-1));
		int cat1 = (int)(rnd::GetRandom().Uniform() * (GetNcomponent()-2));  
		int cat2 = cat1 + 1;
		double logMetropolis = (occupancy[cat1] * log(1 - V[cat2])) - (occupancy[cat2] * log(1-V[cat1]));
		int accepted = (log(rnd::GetRandom().Uniform()) < logMetropolis);
		if (accepted)	{
			total += 1.0;
			SwapComponents(cat1,cat2);
		}
	}

	return total /= nrep;
}
*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
