
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#include "StringStreamUtils.h"
#include <cassert>
#include "RASCATGammaPhyloProcessVI.h"
#include "Parallel.h"
#include <string>

void RASCATGammaPhyloProcessVI::GlobalUpdateParameters()	{
		// MPI2
	// should send the slaves the relevant information
	// about model parameters

	// for this model, should broadcast
	// double alpha
	// int Ncomponent
	// int* alloc
	// double* rr
	// double** profile
	// double* brancharray
	// (but should first call PutBranchLengthsIntoArray())
	// 
	// upon receiving this information
	// slave should 
	// store it in the local copies of the variables
	// and then call
	// SetBranchLengthsFromArray()
	// SetAlpha(inalpha)

	assert(myid == 0);

	// ResampleWeights();
	RenormalizeProfiles();

	int i,j,nbranch = GetNbranch(),ncat = GetNsite(),ni,nd,L1,L2;
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 3 + nbranch + ncat + L1*L2 + L1*L2 + 1;
	ni = 1 + GetNsite();
	int ivector[ni];
	double dvector[nd];
	MESSAGE signal = PARAMETER_DIFFUSION;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	// First we assemble the vector of doubles for distribution
	int index = 0;
	dvector[index] = alpha;
	index++;
	dvector[index] = branchalpha;
	index++;
	dvector[index] = branchbeta;
	index++;
        /*for(i=0; i<ncat; ++i)  {
                dvector[index] = RateAlphaVI[i];
                index++;
        }
        for(i=0; i<ncat; ++i)  {
                dvector[index] = RateBetaVI[i];
                index++;
        }
        for(i=0; i<nbranch; ++i) {
		dvector[index] = branchalphaVI[i];
		index++;
	}
        for(i=0; i<nbranch; ++i) {
		dvector[index] = branchbetaVI[i];
		index++;
	}*/
	for(i=0; i<nbranch; ++i) {
		dvector[index] = blarray[i];
		index++;
	}
	for(i=0; i<ncat; ++i)  {
                dvector[index] = rate[i];
        }
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = profile[i][j];
			index++;
		}
	}
        for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dvector[index] = dirweightVI[i][j];
			index++;
		}
	}
	dvector[index] = kappa;
	index++;
        /*for(i=0; i<L1; ++i)  {
                dvector[index] = kappa_alpha[i];
		index++;
        }
        for(i=0; i<L1; ++i)  {
                dvector[index] = kappa_beta[i];
		index++;
        }*/
	// Now the vector of ints
	ivector[0] = GetNcomponent();
	for(i=0; i<GetNsite(); ++i) {
		ivector[1+i] = DPProfileProcessVI::alloc[i];
	}

	// Now send out the doubles and ints over the wire...
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void RASCATGammaPhyloProcessVI::SlaveExecute(MESSAGE signal)	{

	assert(myid > 0);

	switch(signal) {

	/*
	case PRINT_TREE:
		SlavePrintTree();
		break;
	*/
        case UPDATE_ESTBRANCH:
		SlaveEstimateBranchParameter();
		break;
        case UPDATE_MOVEBRANCH:
		SlaveMoveBranchParamter();
		break;
	case UPDATE_ESTRATE:
		SlaveEstimateRateParameter();
		break;
        case UPDATE_MOVERATE:
		SlaveMoveRateParameter();
		break;
        case UPDATE_MOVEBRAN:
                SlaveMoveLength();
                break;
        case UPDATE_MOVEPRRATE:
                SlaveMoveRate();
                break;
        case MEANRATEVI:
                SlaveGetMeanRateVI();
                break;
        case ELBORATEVI:
                SlaveELBORate();
                break;
        case RATEALPHAVI:
                SlaveGetRateAlphaVI();
                break;
        case RATEBETAVI:
                SlaveGetRateBetaVI();
                break;
	default:
		PhyloProcess::SlaveExecute(signal);
	}
}


void RASCATGammaPhyloProcessVI::SlaveUpdateParameters()	{
	int i,j,L1,L2,ni,nd,nbranch = GetNbranch(),ncat = GetNsite();
	L1 = GetNmodeMax();
	L2 = GetDim();
	nd = 3 + nbranch + ncat + L1*L2 + L1*L2 + 1;
	ni = 1 + GetNsite();
	int* ivector = new int[ni];
	double* dvector = new double[nd];
	MPI_Bcast(ivector,ni,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dvector,nd,MPI_DOUBLE,0,MPI_COMM_WORLD);

	int index = 0;
        alpha = dvector[index];
        index++;
        branchalpha = dvector[index];
        index++;
        branchbeta = dvector[index];
        index++;
	/*for(i=0; i<ncat; ++i)  {
                RateAlphaVI[i] = dvector[index];
                index++;
        }
        for(i=0; i<ncat; ++i)  {
                RateBetaVI[i] = dvector[index];
                index++;
        }
        for(i=0; i<nbranch; ++i) {
		branchalphaVI[i] = dvector[index];
		index++;
	}
        for(i=0; i<nbranch; ++i) {
		branchbetaVI[i] = dvector[index];
		index++;
	}*/
	for(i=0; i<nbranch; ++i) {
		blarray[i] = dvector[index];
		index++;
	}
        for(i=0; i<ncat; ++i)  {
                rate[i] = dvector[index];
        }
	for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			profile[i][j] = dvector[index];
			index++;
		}
	}
        for(i=0; i<L1; ++i) {
		for(j=0; j<L2; ++j) {
			dirweightVI[i][j] = dvector[index];
			index++;
		}
	}
        kappa = dvector[index];
        index++;
	/*for(i=0; i<L1; ++i)  {
                kappa_alpha[i] = dvector[index];
		index++;
        }
        for(i=0; i<L1; ++i)  {
                kappa_beta[i] = dvector[index];
		index++;
        }*/

	Ncomponent = ivector[0];
	for(i=0; i<GetNsite(); ++i) {
		DPProfileProcessVI::alloc[i] = ivector[1+i];
	}
	delete[] dvector;
	delete[] ivector;

	UpdateZip();
}

/*double RASCATGammaPhyloProcessVI::CheckELBO()  {

       GetTotalSitePCAT();
       
       ELBO += ELBORate() + ELBOLength() + ELBOWeight() + ELBOkappa() + ELBODirWeight(); 
       for (int k=0; k<GetNcomponent(); k++)  {
            for (int i=0; i<GetDim(); i++)       {
                 ELBO += (GetPCAT(k) * (gsl_sf_psi(dirweightVI[k][i]) - TotDirWeightVI(k,i)));
            }
       }
       for (int i=GetSiteMin(); i<GetSiteMax(); i++)  {
            ELBO += (GetSiteRateSuffStatCount(i) * (gsl_sf_psi(GetRateAlphaVI(i)) - log(GetRateBetaVI(i))));     
       }
       for (int j=1; j<GetNbranch(); j++)  {
            ELBO += (GetBranchLengthSuffStatCount(j) * (gsl_sf_psi(GetbranchalphaVI(j)) - log(GetbranchbetaVI(j))));
       }
       for (int i=GetSiteMin(); i<GetSiteMax(); i++)   {
            for (int j=1; j<GetNbranch(); j++)              {
                 ELBO -= ((GetbranchalphaVI(j) / GetbranchbetaVI(j)) * (GetRateAlphaVI(i) / GetRateBetaVI(i)));
            }
       }
       return ELBO;
}*/

void RASCATGammaPhyloProcessVI::ReadPB(int argc, char* argv[])	{

	string name = "";

	int burnin = -1;
	int every = 1;
	int until = -1;
	int ppred = 0;
	int ss = 0;
	int map = 0;
	// 1 : plain ppred (outputs simulated data)
	// 2 : diversity statistic
	// 3 : compositional statistic
	int cv = 0;
	int sitelogl = 0;
	int rates = 0;
	string testdatafile = "";
	int rateprior = 0;
	int profileprior = 0;
	int rootprior = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-div")	{
				ppred = 2;
			}
			else if (s == "-comp")	{
				ppred = 3;
			}
			else if (s == "-ppred")	{
				ppred = 1;
			}
			else if (s == "-ppredrate")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rateprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rateprior = 0;
				}
				else	{
					cerr << "error after ppredrate: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredprofile")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					profileprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					profileprior = 0;
				}
				else	{
					cerr << "error after ppredprofile: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-ppredroot")	{
				i++;
				string tmp = argv[i];
				if (tmp == "prior")	{
					rootprior = 1;
				}
				else if ((tmp == "posterior") || (tmp == "post"))	{
					rootprior = 0;
				}
				else	{
					cerr << "error after ppredroot: should be prior or posterior\n";
					throw(0);
				}
			}
			else if (s == "-sitelogl")	{
				sitelogl = 1;
			}
			else if (s == "-r")	{
				rates = 1;
			}
			else if (s == "-cv")	{
				cv = 1;
				i++;
				testdatafile = argv[i];
			}
			else if (s == "-ss")	{
				ss = 1;
			}
			else if (s == "-map")	{
				map = 1;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					string tmp = argv[i];
					if (IsInt(tmp))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else {
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	if (until == -1)	{
		until = GetSize();
	}
	if (burnin == -1)	{
		burnin = GetSize() / 5;
	}

	if ((GetNprocs() == 1) && (ppred || cv || sitelogl))	{
		cerr << "error : should run readpb_mpi in mpi mode, with at least 2 processes\n";
		MPI_Finalize();
		exit(1);
	}

	if (ss)	{
		ReadSiteProfiles(name,burnin,every,until);
	}
	else if (map)	{
		ReadMap(name,burnin,every,until);
	}
	else if (cv)	{
		ReadCV(testdatafile,name,burnin,every,until);
	}
	else if (sitelogl)	{
		ReadSiteLogL(name,burnin,every,until);
	}
	else if (rates)	{
		ReadSiteRates(name,burnin,every,until);
	}
	else if (ppred)	{
		PostPred(ppred,name,burnin,every,until,rateprior,profileprior,rootprior);
	}
	else	{
		Read(name,burnin,every,until);
	}
}

void RASCATGammaPhyloProcessVI::ReadSiteProfiles(string name, int burnin, int every, int until)	{

	double** sitestat = new double*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		sitestat[i] = new double[GetDim()];
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] = 0;
		}
	}
	ifstream is((name + ".chain").c_str());
	if (!is)	{
		cerr << "error: no .chain file found\n";
		exit(1);
	}

	cerr << "burnin : " << burnin << "\n";
	cerr << "until : " << until << '\n';
	int i=0;
	while ((i < until) && (i < burnin))	{
		FromStream(is);
		i++;
	}
	int samplesize = 0;

	while (i < until)	{
		cerr << ".";
		cerr.flush();
		samplesize++;
		FromStream(is);
		i++;

		for (int i=0; i<GetNsite(); i++)	{
			double* p = GetProfile(i);
			for (int k=0; k<GetDim(); k++)	{
				sitestat[i][k] += p[k];
			}
		}
		int nrep = 1;
		while ((i<until) && (nrep < every))	{
			FromStream(is);
			i++;
			nrep++;
		}
	}
	cerr << '\n';
	
	ofstream os((name + ".siteprofiles").c_str());
	for (int k=0; k<GetDim(); k++)	{
		os << GetStateSpace()->GetState(k) << ' ';
	}
	os << '\n';
	os << '\n';
	for (int i=0; i<GetNsite(); i++)	{
		os << i + 1;
		for (int k=0; k<GetDim(); k++)	{
			sitestat[i][k] /= samplesize;
			os << '\t' << sitestat[i][k];
		}
		os << '\n';
	}

	cerr << "mean site-specific profiles in " << name << ".siteprofiles\n";
	cerr << '\n';
}

