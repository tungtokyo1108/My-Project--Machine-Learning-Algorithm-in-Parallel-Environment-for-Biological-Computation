
/********************


**********************/


#include "PoissonSubstitutionProcess.h"

#include "Parallel.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Poisson Substitution Processes
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// same methods as above, but specialized for the CAT poisson model
// probably less interesting to parallelize, at least in a first step


//-------------------------------------------------------------------------
//	* conditional likelihood propagation
//	(CPU level 3)
//-------------------------------------------------------------------------

/*void PoissonSubstitutionProcess::Propagate(double*** from, double*** to, double time, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		const double* stat = GetStationary(i);
		for (int j=0; j<GetNrate(i); j++)	{
			if ((! condalloc) || (ratealloc[i] == j))	{
				double* tmpfrom = from[i][j];
				double* tmpto = to[i][j];
				double expo = exp(-GetRate(i,j) * time);
				double tot = 0;
				int nstate = GetNstate(i);
				for (int k=0; k<nstate; k++)	{
					tot += (*tmpfrom++) * (*stat++);
					// tot += tmpfrom[k] * stat[k];
				}
				tmpfrom -= nstate;
				stat -= nstate;
				tot *= (1-expo);
				for (int k=0; k<nstate; k++)	{	
					(*tmpto++) = expo * (*tmpfrom++) + tot;
					// tmpto[k] = expo * tmpfrom[k] + tot;
				}
				(*tmpto) = (*tmpfrom);
				tmpto -= nstate;
				tmpfrom -= nstate;
				// tmpto[GetNstate(i)] = tmpfrom[GetNstate(i)];
			}
		}
	}
}*/

/*void PoissonSubstitutionProcess::Propagate(double** from, double** to, double time, bool condalloc)	{
	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		double* tmpfrom = from[i];
		double* tmpto = to[i];
		double expo = exp(-GetRate(i) * time);
		double tot = 0;
		int nstate = GetNstate(i);
		for (int k=0; k<nstate; k++)	{
			tot += (*tmpfrom++) * (*stat++);
		}
		tmpfrom -= nstate;
		stat -= nstate;
		tot *= (1-expo);
		for (int k=0; k<nstate; k++)	{	
			(*tmpto++) = expo * (*tmpfrom++) + tot;
		}
		(*tmpto) = (*tmpfrom);
		tmpto -= nstate;
		tmpfrom -= nstate;
	}
}*/

void PoissonSubstitutionProcess::Propagate(double** from, double** to, double time, bool condallc)  {

       for (int i=sitemin; i<sitemax; i++)  {
               const double* stat = GetStationary(i);
               double* tmpfrom = from[i];
               double* tmpto = to[i];
               double expo = exp(-GetRate(i) * time);
               double tot = 0;
               int nstate = GetNstate(i);
               /*for (int k=0; k<nstate; k++)   {
                       tot += (*stat++);
               }
               stat -= nstate;
               tot *= (1-expo);
               for (int k=0; k<nstate; k++)   {
                       (*tmpto++) = tot;
               }
               tmpto -= nstate;*/

               double remanence = (1-expo) * (*stat++);
               for (int k=0; k<nstate; k++)   {
                       (*tmpto++) *= remanence; 
               }
               tmpto -= nstate;
       }
}

/*
// version directly unzipped
void PoissonSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{

	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetProfile(i);
		int j = ratealloc[i];
		double expo = exp(-GetRate(i,j) * time);
		if (rnd::GetRandom().Uniform() < expo)	{
			statedown[i] = stateup[i];
		}
		else	{
			double u = rnd::GetRandom().Uniform();
			int k = 0;
			double cumul = stat[k];
			while ((k < GetDim()) && (u > cumul))	{
				k++;
				if (k == GetDim())	{
					cerr << "error in PoissonSubstitutionProcess::SimuPropagate: overflow\n";
					exit(1);
				}
				cumul += stat[k];
			}
			statedown[i] = k;
		}
	}
}
*/

/*void PoissonSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{

	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		int nstate = GetNstate(i);
		int j = ratealloc[i];
		double expo = exp(-GetRate(i,j) * time);
		if (rnd::GetRandom().Uniform() < expo)	{
			statedown[i] = stateup[i];
		}
		else	{
			double u = rnd::GetRandom().Uniform();
			int k = 0;
			double cumul = stat[k];
			while ((k < nstate) && (u > cumul))	{
				k++;
				if (k == GetDim())	{
					cerr << "error in PoissonSubstitutionProcess::SimuPropagate: overflow\n";
					exit(1);
				}
				cumul += stat[k];
			}
			statedown[i] = k;
		}
	}
}*/

void PoissonSubstitutionProcess::SimuPropagate(int* stateup, int* statedown, double time)	{
       for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		int nstate = GetNstate(i);
		double expo = exp(-GetRate(i) * time);
		if (rnd::GetRandom().Uniform() < expo)	{
			statedown[i] = stateup[i];
		}
		else	{
			double u = rnd::GetRandom().Uniform();
			int k = 0;
			double cumul = stat[k];
			while ((k < nstate) && (u > cumul))	{
				k++;
				if (k == GetDim())	{
					cerr << "error in PoissonSubstitutionProcess::SimuPropagate: overflow\n";
					exit(1);
				}
				cumul += stat[k];
			}
			statedown[i] = k;
		}
	}
}

//-------------------------------------------------------------------------
//	* sample substitution mappings conditional on states at nodes 
//	(CPU level 3)
//-------------------------------------------------------------------------

// root version
BranchSitePath** PoissonSubstitutionProcess::SampleRootPaths(int* state)	{
	// BranchSitePath** patharray = new BranchSitePath*[sitemax - sitemin];
	BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		patharray[i] = new BranchSitePath(0,state[i]);
	}
	return patharray;
}

// general version
/*BranchSitePath** PoissonSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time) 	{
	BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
	for (int i=sitemin; i<sitemax; i++)	{
		const double* stat = GetStationary(i);
		double rate = GetRate(i);
		double l = rate * time;
		int dup = stateup[i];
		int ddown = statedown[i];
		double pi = stat[ddown];
                double clusterprofile = GetPCAT(ddown);
		int m = 0;
		int mmax = 1000;
		
		if (dup == ddown)	{
			double fact = pi * exp(-l);
			double total = exp(-l);
			double q = rnd::GetRandom().Uniform() * (exp(-l) * (1 - pi) + pi);
			while ((m<mmax) && (total < q))	{
				m++;
				fact *= l / m;
				total += fact;
			}
			if (m == mmax)	{
				suboverflowcount ++;
			}
		}
		else	{
			double fact = pi * exp(-l);
			double total = 0;
			double q = rnd::GetRandom().Uniform() * (1 - exp(-l)) * pi;
			while ((m<mmax) && (total < q))	{
				m++;
				fact *= l / m;
				total += fact;
			}
			if (m == mmax)	{
				suboverflowcount ++;
			}
		}
		patharray[i] = new BranchSitePath(m,ddown);
	}
	return patharray;
}*/

BranchSitePath** PoissonSubstitutionProcess::SamplePathsVI(int* stateup, int* statedown, double timealpha, double timebeta)   {

        BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
        for (int i=sitemin; i<sitemax; i++)   {
             const double* stat = GetStationaryDirweightVI(i);
             double ratealpha = GetRateAlphaVI(i);
             double ratebeta = GetRateBetaVI(i);
             int dup = stateup[i];
             int ddown = statedown[i];
             double pi = stat[ddown];
             double pidirweight = 0;
             for (int k=0; k<GetDim(); k++)     {
                  pidirweight += stat[k];
             }
             int m = 0;
             if (dup == ddown)   {
                    m += gsl_sf_psi(pi) - gsl_sf_psi(pidirweight) + gsl_sf_psi(ratealpha) - log(ratebeta) + gsl_sf_psi(timealpha) - log(timebeta) - ((ratealpha / ratebeta)*(timealpha / timebeta));
                    m -= log(exp(-1/10) + (1 - exp(-1/10) * 1/20)) + ((1/20 -1) / (1/20 * exp(-1/10) - 1/20 + 1)) * ((ratealpha / ratebeta) * (timealpha / timebeta) - 1/10);
                    m -= ((exp(-1/10) - 1) / ((exp(-1/10) - 1) * 1/20 + 1)) * ((pi / pidirweight) - 1/20);
             }
             else   {
                    m = 1;
             }
             patharray[i] = new BranchSitePath(m,ddown);
        }
        return patharray;
}

BranchSitePath** PoissonSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time)   {

       BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
       for (int i=sitemin; i<sitemax; i++)   {
              const double* stat = GetStationary(i);
              // const double* stat = GetStationaryDirweightVI(i);
              double rate = GetRate(i);
              double l = rate * time;
              int dup = stateup[i];
              int ddown = statedown[i];
              double pidown = stat[ddown];
              double piup = stat[dup];
              double norm = 0;
              double total = 0;
 
              if (dup == ddown)      {
                     double w = (exp(-l) + (1 - exp(-l)) * piup) * pidown;
                     norm += w;
                     total += w * l * pidown / (exp(-l) + (1 - exp(-l)) * pidown);                    
              }
              else   {
                     double w = (1 - exp(-l)) * pidown * piup;
                     norm += w;
                     total += w * l / (1 - exp(-l));
              }
              int m = total / norm;
              patharray[i] = new BranchSitePath(m,ddown);
       }
       return patharray;
}

/*BranchSitePath** PoissonSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time)  {

        BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
        for (int i=sitemin; i<sitemax; i++)   {
             for (int j=1; j<GetNbranch(); j++)   {
                     const double* stat = GetStationaryDirweightVI(i);
                     int dup = stateup[i];
                     int ddown = statedown[i];
                     double pidown = stat[ddown];
                     double piup = stat[dup];
                     double norm = 0;
                     double total = 0;

                     if (dup == ddown)   {
                           double w = (exp(-(GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j))) + (1 - exp(-(GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)))) * piup) * pidown;
                           norm += w;
                           total += w * (GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)) * pidown / (exp(-(GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j))) + (1 - exp(-(GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)))) * pidown);
                     }
                     else  {
                           double w = (1 - exp(- (GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)))) * pidown * piup;
                           norm += w;
                           total += w * (GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)) / (1 - exp(- (GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j))));
                     }
                     int m = total / norm;
                     patharray[i] = new BranchSitePath(m,ddown);     
             }  
        }
        return patharray;
}*/

/*BranchSitePath** PoissonSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time)   {

        BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
        for (int i=sitemin; i<sitemax; i++)    {
             const double* stat = GetStationaryDirweightVI(i);
             int dup = stateup[i];
             int ddown = statedown[i];
             double pidown = stat[ddown];
             double piup = stat[dup];
             int m = 0;
             for (int j=1; j<GetNbranch(); j++)    {    
                  double pidirweight = 0;
                  for (int k=0; k<GetDim(); k++)  {
                       pidirweight += stat[k];
                  }
                  if (dup == ddown)     {
                        // m += gsl_sf_psi(pidown) - gsl_sf_psi(pidirweight) + gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i)) + gsl_sf_psi(GetbranchalphaVI(j)) - log(GetbranchbetaVI(j)) - ((GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)));
                        m += log(exp(-1/10) + (1 - exp(-1/10) * 1/20)) + ((1/20 -1) / (1/20 * exp(-1/10) - 1/20 + 1)) * ((GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)) - 1/10);
                        m += ((exp(-1/10) - 1) / ((exp(-1/10) - 1) * 1/20 + 1)) * ((pidown / pidirweight) - 1/20);
                  }
                  else    {
                        m = 1;
                  }
             }
             patharray[i] = new BranchSitePath(m,ddown);
        }
        return patharray;
}*/

/*BranchSitePath** PoissonSubstitutionProcess::SamplePaths(int* stateup, int* statedown, double time)    {

        BranchSitePath** patharray = new BranchSitePath*[GetNsite()];
        for (int i=sitemin; i<sitemax; i++)    {
             const double* stat = GetStationaryDirweightVI(i);
             int dup = stateup[i];
             int ddown = statedown[i];
             double pidown = stat[ddown];
             double piup = stat[dup];
             int m = 0;
             double integral = 0;
             int N = 10000;
             for (int j=1; j<GetNbranch(); j++)   {
                  double pidirweight = 0;
                  for (int k=0; k<GetDim(); k++)     {
                       pidirweight += stat[k];
                  }
                  if (dup == ddown)     {
                       m += gsl_sf_psi(pidown) - gsl_sf_psi(pidirweight) + gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i)) + gsl_sf_psi(GetbranchalphaVI(j)) - log(GetbranchbetaVI(j)) - ((GetRateAlphaVI(i) / GetRateBetaVI(i)) * (GetbranchalphaVI(j) / GetbranchbetaVI(j)));
                       for (int i=0; i<N; i++)    {
                            for (int j=0; j<N; j++)   {
                                 for (int k=0; k<N; k++)    {
                                       integral = exp((pidown - 1) * log(i/N) + rnd::GetRandom().logGamma(pidirweight) - rnd::GetRandom().logGamma(pidown)) * exp( GetRateAlphaVI(i) * log(GetRateBetaVI(i)) - rnd::GetRandom().logGamma(GetRateAlphaVI(i)) + (GetRateAlphaVI(i) - 1) * log(j/N) - GetRateBetaVI(i) * (j/N)) * exp( GetbranchalphaVI(j) * log(GetbranchbetaVI(j)) - rnd::GetRandom().logGamma(GetbranchalphaVI(j)) + (GetbranchalphaVI(j) - 1) * log(k/N) - GetbranchbetaVI(j) * (k/N)) * log( exp(-j*k/N^2) + (1 - exp(-j*k/N^2)) * i/N );
                                 }
                            }
                       }
                       integral /= N^3;
                       m += integral;
                  }
                  else    {
                       m = 1;
                  }
             }
             patharray[i] = new BranchSitePath(m,ddown);
        }
        return patharray;
}*/

//-------------------------------------------------------------------------
//	* gather sufficient statistics 
//	(CPU level 3)
//-------------------------------------------------------------------------


void PoissonSubstitutionProcess::AddSiteRateSuffStat(int* siteratesuffstatcount, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		siteratesuffstatcount[i] += patharray[i]->GetNsub();
	}
}

void PoissonSubstitutionProcess::AddBranchLengthSuffStat(int& count, BranchSitePath** patharray)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		count += patharray[i]->GetNsub();
	}
}

void PoissonSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root)	{
	for (int i=sitemin; i<sitemax; i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		if (root || patharray[i]->GetNsub())	{
			siteprofilesuffstatcount[i][GetRandomStateFromZip(i,patharray[i]->GetFinalState())]++;
		}
	}
}

/*void PoissonSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, bool root)  {

        for (int i=sitemin; i<sitemax; i++)   {
                if (root || patharray[i]->GetNsub())   {
                        siteprofilesuffstatcount[i][GetRandomStateFromZip(i,patharray[i]->GetFinalState())];
                }
        }
}*/

/*void PoissonSubstitutionProcess::AddSiteProfileSuffStat(int** siteprofilesuffstatcount, BranchSitePath** patharray, double numbertypesubstitution, bool root)    {

       for (int i=sitemin; i<sitemax; i++)    {
              if (root || patharray[i]->GetNsub())   {
                    for (int k=0; k<GetDim(); k++)     {
                             siteprofilesuffstatcount[i][k] = numbertypesubstitution;
                    }
              }
       }
}*/

void PoissonSubstitutionProcess::ChooseTrueStates(BranchSitePath** patharray, int* nodestateup, int* nodestatedown, bool root)	{
	for (int i=sitemin; i<sitemax; i++)	{
		int tmp = nodestateup[i];
		if (root || patharray[i]->GetNsub())	{
			tmp = GetRandomStateFromZip(i,patharray[i]->GetFinalState());
		}
		//cerr  << ' '<< i  << ' '<< patharray[i]->GetNsub() << ' '<< tmp  << '\n';
		nodestatedown[i] = tmp;
	}
}

void PoissonSubstitutionProcess::ChooseRootTrueStates(int* nodestate)	{

	for (int i=sitemin; i<sitemax; i++)	{
		int tmp = GetRandomStateFromZip(i,nodestate[i]);
		nodestate[i] = tmp;
	}
}

//-------------------------------------------------------------------------
//	* recomputing the equilibrium frequency profiles of the recoded process
//
//-------------------------------------------------------------------------


void PoissonSubstitutionProcess::CreateZip()	{
	if (! zipstat)	{
		zipstat = new double*[GetNsite()];
		for (int i=0; i<GetNsite(); i++)	{
			zipstat[i] = new double[GetDim()];
			// zipstat[i] = new double[GetNstate(i)];
		}
	}
}

void PoissonSubstitutionProcess::DeleteZip()	{
	for (int i=0; i<GetNsite(); i++)	{
		delete[] zipstat[i];
	}
	delete[] zipstat;
	zipstat = 0;
}

void PoissonSubstitutionProcess::UpdateZip()	{
	
	for (int i=0; i<GetNsite(); i++)	{
		UpdateZip(i);
	}
}

void PoissonSubstitutionProcess::UpdateZip(int i)	{
		double total = 0;
		double* pi = GetProfile(i);
		for (int k=0; k<GetOrbitSize(i); k++)	{
			int n = GetStateFromZip(i,k);
			zipstat[i][k] = 0;
			zipstat[i][k] = pi[GetStateFromZip(i,k)];
			total += zipstat[i][k];
		}
		if (GetZipSize(i) > GetOrbitSize(i))	{
			zipstat[i][GetOrbitSize(i)] = 1-total;
		}
}

/*int PoissonSubstitutionProcess::GetRandomStateFromZip(int site, int zipstate)	{
	int truestate = 0;
	if ((GetZipSize(site) != GetOrbitSize(site)) && (zipstate == GetOrbitSize(site)))	{
		double v = rnd::GetRandom().Uniform();
		double u = zipstat[site][GetOrbitSize(site)] * v;
		double total = 0;
		double* pi = GetProfile(site);
		int choose = -1;
		while ((choose < GetDim()) && (total < u))	{
			choose ++;
			if (choose == GetDim())	{
				cerr << "error in getstatefromzip\n";
				cerr << choose << '\t' << GetDim() << '\n';
				cerr << total << '\n';
				cerr << v << '\t' << zipstat[site][GetOrbitSize(site)] << '\t' << u << '\n';
				cerr << total - zipstat[site][GetOrbitSize(site)] << '\n';
				double newtotal = 0;
				cerr << '\n';
				for (int k=0; k<GetDim(); k++)	{
					cerr << pi[k] << '\t';
					cerr << InOrbit(site,k) << '\n';
					if (! InOrbit(site,k))	{
						newtotal += pi[k];
					}
				}
				cerr << '\n';
				for (int k=0; k<=GetOrbitSize(site); k++)	{
					cerr << zipstat[site][k] << '\n';
				}
				cerr << "new total : " << newtotal << '\t' << newtotal - total << '\t' << total - choose << '\n';
				exit(1);
			}
			if (!InOrbit(site,choose))	{
				total += pi[choose];
			}
		}
		truestate = choose;
	}
	else	{
		truestate = GetStateFromZip(site,zipstate);
	}
	return truestate;
}*/

int PoissonSubstitutionProcess::GetRandomStateFromZip(int site, int zipstate)   {
 
       int truestate = 0;
       if ((GetZipSize(site) != GetOrbitSize(site)) && (zipstate == GetOrbitSize(site)))   {
              double* pi = GetStationaryDirweightVI(site);
              double ratealpha = GetRateAlphaVI(site);
              double ratebeta = GetRateBetaVI(site);
              double total = 0;
              double norm = 0;
              for (int j=1; j<GetNbranch(); j++)   {
                   double branchalpha = GetbranchalphaVI(j);
                   double branchbeta = GetbranchbetaVI(j);       
                   for (int k=0; k<GetDim(); k++)      {
                      double w = exp(-(ratealpha / ratebeta) * (branchalpha / branchbeta)) + (1 - exp(-(ratealpha / ratebeta) * (branchalpha / branchbeta))) * pi[k];
                      norm += w;
                      total += w * (1 - exp(-(ratealpha / ratebeta) * (branchalpha / branchbeta))) * pi[k] / (exp(-(ratealpha / ratebeta) * (branchalpha / branchbeta)) + (1 - exp(-(ratealpha / ratebeta) * (branchalpha / branchbeta))) * pi[k]);
                   }
              }
              truestate = total / norm;
       }
       return truestate;
}

void PoissonSubstitutionProcess::UnzipBranchSitePath(BranchSitePath** patharray, int* nodestateup, int* nodestatedown){
	for (int i=sitemin; i<sitemax; i++)	{
		int nsub = patharray[i]->GetNsub();
		patharray[i]->nsub=0;
		double* times = new double[nsub+1];
		for(int j = 0; j < nsub; j++){
			times[j] = rnd::GetRandom().Uniform();
			for(int k = 0; k < j; k++){
				if(times[k]>times[j]	){
					times[nsub] = times[k];
					times[k] = times[j];
					times[j] = times[nsub];
				}
			}
		}
		double mem = 0;
		for(int j = 0; j < nsub; j++){
			times[nsub] = times[j];
			times[j] = times[j] - mem;
			mem = times[nsub];
		}
		times[nsub]=1-mem;

		int previousstate = nodestateup[i];
		patharray[i]->Init()->SetState(previousstate);
		double* pi = GetProfile(i);
		for(int j = 0; j < nsub-1; j++){
			int newstate = rnd::GetRandom().DrawFromDiscreteDistribution(pi, GetDim());
			if(newstate != previousstate){
			      patharray[i]->Append(newstate, times[j]);
			      previousstate=newstate;
			}
			else{
			      times[j+1]+=times[j];
			}
		}
		if(previousstate == nodestatedown[i]){
			times[nsub] += times[nsub-1];
		}
		else{
			patharray[i]->Append(nodestatedown[i], times[nsub-1]);
		}
		patharray[i]->Last()->SetRelativeTime(times[nsub]);
		delete[] times;
	}
}

