
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#include "PoissonSBDPProfileProcessVI.h"
#include <cassert>
#include "Parallel.h"

void PoissonSBDPProfileProcessVI::Create(int innsite, int indim)  {
              
                PoissonDPProfileProcessVI::Create(innsite,indim);
		SBDPProfileProcessVI::Create(innsite,indim);
                UpdateProCAT = new double[GetNmodeMax()];
                /*PCAT = new double*[GetNsite()];
                for (int site = 0; site < GetNsite(); site++) {
                       PCAT[site] = new double[GetNmodeMax()];
                }
                TransPCAT = new double*[GetNmodeMax()];
                for (int j = 0; j < GetNmodeMax(); j++)  {
                       TransPCAT[j] = new double[GetNsite()];
                }*/
                NormalizePCAT = new double*[GetNsite()];
                for (int site = 0; site < GetNsite(); site++) {
                       NormalizePCAT[site] = new double[GetNmodeMax()];
                }

                /*for (int mode = 0; mode < GetNmodeMax(); mode++) {
                     
                     // PCAT[mode] = rnd::GetRandom().Uniform();
                     PCAT[mode] = 0.5;
                     // PCAT[mode] = exp(GetProfileParameter(mode));
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

double PoissonSBDPProfileProcessVI::GlobalMixMove(int nrep, int nallocrep, double epsilon)	{

	if (Ncomponent != GetNmodeMax())	{
		cerr << "error in sbdp inc dp move: number of components\n";
		exit(1);
	}
         
	int K0 = GetNmodeMax();

	// send mixmove signal and tuning parameters
	MESSAGE signal = MIX_MOVE;
	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
	int itmp[3];
	itmp[0] = nrep;
	itmp[1] = nallocrep;
	itmp[2] = K0;
	MPI_Bcast(itmp,3,MPI_INT,0,MPI_COMM_WORLD);

	// split Nsite among GetNprocs()-1 slaves
	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
	}

	double* tmp = new double[Ncomponent * GetDim() + 1];	
        double* tmpVI = new double[GetNsite() * GetNmodeMax() + 1]; 
        for (int rep=0; rep<nrep; rep++) {
                
                ResampleWeights();
                MPI_Bcast(weightVI,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// here slaves do realloc moves

		// mpi receive new allocations
		MPI_Status stat;
                // double tmpproVI[GetNsite() * GetNmodeMax()+1];
		int tmpalloc[GetNsite()+1];
		for(int i=1; i<GetNprocs(); ++i) { 
                        MPI_Recv(tmpVI,(smax[i-1]-smin[i-1])*GetNmodeMax()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);     
			MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
                        int t = 0;
			for(int j=smin[i-1]; j<smax[i-1]; ++j) {
				alloc[j] = tmpalloc[j];
				if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))	{
					cerr << "alloc overflow\n";
					exit(1);
				}
                                
                                for(int i=0; i<GetNmodeMax(); i++) {
                                    NormalizePCAT[j][i] = tmpVI[t];
                                    t++;
                                }
			}
		}

		// MPI_Barrier(MPI_COMM_WORLD);

		// broadcast new allocations
                // MPI_Bcast(NormalizePCAT,GetNsite() * GetNmodeMax()+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// here slaves do profile moves

		UpdateOccupancyNumbers();

		// collect final values of profiles (+ total acceptance rate) from slaves

		// split Nmode among GetNprocs()-1 slaves
		int mwidth = GetNcomponent()/(GetNprocs()-1);
		int mmin[GetNprocs()-1];
		int mmax[GetNprocs()-1];
		for(int i=0; i<GetNprocs()-1; ++i) {
			mmin[i] = mwidth*i;
			mmax[i] = mwidth*(1+i);
			if (i == (GetNprocs()-2)) mmax[i] = GetNcomponent();
		}

		MPI_Status stat2;
		double total = 0;
		for(int i=1; i<GetNprocs(); ++i) {
			MPI_Recv(tmp,(mmax[i-1]-mmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat2);
			int l = 0;
			for(int j=mmin[i-1]; j<mmax[i-1]; ++j) {
				for (int k=0; k<GetDim(); k++)	{
					profile[j][k] = tmp[l];
					l++;
				}
			}
			total += tmp[l]; // (sum all acceptance rates)
		}

		// MPI_Barrier(MPI_COMM_WORLD);
		// resend all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        
        delete[] tmpVI;
	delete[] tmp;

	return 1;
}

/*double PoissonSBDPProfileProcessVI::GlobalMixMove(int nrep, int nallocrep, double epsilon)   {

       if (Ncomponent != GetNmodeMax())  {
               cerr << "error in sbdp in dp move: number of components\n";
               exit(1);
       }

       int K0 = GetNmodeMax();
       MESSAGE signal = MIX_MOVE;
       MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
       int itmp[3];
       itmp[0] = nrep;
       itmp[1] = nallocrep;
       itmp[2] = K0;
       MPI_Bcast(itmp,3,MPI_INT,0,MPI_COMM_WORLD);

       int width = GetNsite()/(GetNprocs()-1);
       int smin[GetNprocs()-1];
       int smax[GetNprocs()-1];
       for (int i=0; i<GetNprocs()-1; ++i)  {
               smin[i] = width*i;
               smax[i] = width*(1+i);
               if (i == (GetNprocs()-2)) smax[i] = GetNsite();
       }

       double* tmp = new double[Ncomponent * GetDim() +1];
       double* tmpVI = new double[GetNsite() * GetNmodeMax() +1];
       double** TransPCAT = new double*[GetNmodeMax()];
       for (int j=0; j<GetNmodeMax(); j++)  {
               TransPCAT[j] = new double[GetNsite()];
       }
       for (int j=0; j<GetNmodeMax(); j++)  {
               UpdateProCAT[j] = 0;
       }
       for (int rep=0; rep<nrep; rep++)  {
               
               ResampleWeights();
               MPI_Bcast(weightVI,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);
               MPI_Status stat;
               int tmpalloc[GetNsite()+1];
               for (int i=1; i<GetNprocs(); ++i)  {
                        MPI_Recv(tmpVI,(smax[i-1]-smin[i-1])*GetNmodeMax()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
                        MPI_Recv(tmpalloc,GetNsite(),MPI_INT,i,TAG1,MPI_COMM_WORLD,&stat);
                        int t = 0;
                        for (int j=smin[i-1]; j<smax[i-1]; ++j)  {
                                 alloc[j] = tmpalloc[j];
                                 if ((alloc[j] < 0) || (alloc[j] >= Ncomponent))  {
                                         cerr << "alloc overflow\n";
                                         exit(1);
                                 }
                                 
                                 for (int mode=0; mode<GetNmodeMax(); mode++)  {
                                         NormalizePCAT[j][mode] = tmpVI[t];
                                         t++;
                                 }
                        }
                        for (int j=smin[i-1]; j<smax[i-1]; ++j)  {
                                 for (int mode=0; mode<GetNmodeMax(); mode++)   {
                                         TransPCAT[mode][j] = NormalizePCAT[j][mode];
                                 }
                        }
                        for (int mode=0; mode<GetNmodeMax(); mode++)  {
                                 for (int j=smin[i-1]; j<smax[i-1]; j++)  {
                                         UpdateProCAT[mode] += TransPCAT[mode][j];
                                 }
                        }
               }
               MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);
               UpdateOccupancyNumbers();
               
               int mwidth = GetNcomponent()/(GetNprocs()-1);
               int mmin[GetNprocs()-1];
               int mmax[GetNprocs()-1];
               for (int i=0; i<GetNprocs()-1; ++i)  {
                        mmin[i] = mwidth*i;
                        mmax[i] = mwidth*(1+i);
                        if (i == (GetNprocs()-2)) mmax[i] = GetNcomponent();
               }
               
               MPI_Status stat2;
               double total = 0;
               for (int i=1; i<GetNprocs(); ++i)  {
                        MPI_Recv(tmp,(mmax[i-1]-mmin[i-1])*GetDim()+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat2);
                        int l = 0;
                        for (int j=mmin[i-1]; j<mmax[i-1]; ++j)  {
                                for (int k=0; k<GetDim(); k++)     {
                                        profile[j][k] = tmp[l];
                                        l++;
                                }
                        }
                        total += tmp[l];
               }
               MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
       }
       
       delete[] tmpVI;
       delete[] tmp;
       for (int j=0; j<GetNmodeMax(); j++)  {
                delete[] TransPCAT[j];
       }
       delete[] TransPCAT;
    return 1.0;
}*/


void PoissonSBDPProfileProcessVI::SlaveMixMove()	{

	int itmp[3];
	MPI_Bcast(itmp,3,MPI_INT,0,MPI_COMM_WORLD);
	int nrep = itmp[0];
	int nallocrep = itmp[1];
        int K0 = itmp[2];

	// double* mLogSamplingArray = new double[Ncomponent];
        
        PCAT = new double*[GetNsite()];
        for (int site = 0; site < GetNsite(); site++) {
                PCAT[site] = new double[GetNmodeMax()];
        }

	double* tmp = new double[Ncomponent * GetDim() + 1];
        double* tmpVI = new double[GetNsite() * GetNmodeMax() + 1];
	int width = GetNsite()/(GetNprocs()-1);
	int smin[GetNprocs()-1];
	int smax[GetNprocs()-1];
	for(int i=0; i<GetNprocs()-1; ++i) {
		smin[i] = width*i;
		smax[i] = width*(1+i);
		if (i == (GetNprocs()-2)) smax[i] = GetNsite();
	}

        for (int rep=0; rep<nrep; rep++)  {
               
             // ResampleWeights();
             // realloc move
             MPI_Bcast(weightVI,Ncomponent,MPI_DOUBLE,0,MPI_COMM_WORLD);
             // MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);

             /*double totp = 0;
	     for (int mode = 0; mode<K0; mode++)   {
	          totp += weight[mode];
	     }*/

	     /*double totq = 0;
	     for (int mode=K0; mode<GetNmodeMax(); mode++)	{
			totq += weight[mode];
	     }*/
             int t = 0;  
 
	     for (int allocrep=0; allocrep<nallocrep; allocrep++)	{

	          for (int site=smin[GetMyid()-1]; site<smax[GetMyid()-1]; site++)	{                 
                        
                       for (int mode = 0; mode < GetNmodeMax(); mode++)	{
			PCAT[site][mode] = exp(weightVI[mode] + LogStatProb(site,mode));
                        // PCAT[site][mode] = rnd::GetRandom().Uniform();
		        }

                       /*double totalsquare = 0;
                       for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                             totalsquare += (PCAT[site][mode] * PCAT[site][mode]);
                       }
                       double norm = sqrt(totalsquare);*/
                       
                       double total = 0;
                       for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                             total += PCAT[site][mode];
                       }
                       double norm = total;
       
                       for (int mode = 0; mode < GetNmodeMax(); mode++)	{
                              NormalizePCAT[site][mode] = (PCAT[site][mode]) / (norm); 
                              tmpVI[t] = NormalizePCAT[site][mode];
                              t++;  
                       }
                       double MaxProCluster = 0;
                       for (int mode = 0; mode < GetNmodeMax(); mode++)	{
                             if ((!mode) || (MaxProCluster < NormalizePCAT[site][mode]))  {
                               MaxProCluster = NormalizePCAT[site][mode];
                               alloc[site] = mode;
                             }
                       } 
                   }   
              }
                
                MPI_Send(tmpVI,(smax[GetMyid()-1] - smin[GetMyid()-1])*GetNmodeMax() + 1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);            
		MPI_Send(alloc,GetNsite(),MPI_INT,0,TAG1,MPI_COMM_WORLD);

		// MPI_Barrier(MPI_COMM_WORLD);
		// profile move

		// receive new allocations
                // MPI_Bcast(NormalizePCAT,GetNsite() * GetNmodeMax()+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(alloc,GetNsite(),MPI_INT,0,MPI_COMM_WORLD);

		// determine the range of components to move
		UpdateOccupancyNumbers();

		// update sufficient statistics
		UpdateModeProfileSuffStat();

		// split Nmode among GetNprocs()-1 slaves
		int mwidth = GetNcomponent()/(GetNprocs()-1);
		int mmin[GetNprocs()-1];
		int mmax[GetNprocs()-1];
		for(int i=0; i<GetNprocs()-1; ++i) {
			mmin[i] = mwidth*i;
			mmax[i] = mwidth*(1+i);
			if (i == (GetNprocs()-2)) mmax[i] = GetNcomponent();
		}

		double total = 0;
		for (int mode=mmin[GetMyid()-1]; mode<mmax[GetMyid()-1]; mode++)	{
			total += MoveProfile(mode);
		}
		int l = 0;
		for (int mode=mmin[GetMyid()-1]; mode<mmax[GetMyid()-1]; mode++)	{
			for (int k=0; k<GetDim(); k++)	{
				tmp[l] = profile[mode][k];
				l++;
			}
		}
		tmp[l] = total;

		MPI_Send(tmp,(mmax[GetMyid()-1] - mmin[GetMyid()-1])*GetDim()+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
		// MPI_Barrier(MPI_COMM_WORLD);

		// rereceive all profiles
		MPI_Bcast(allocprofile,Ncomponent*GetDim(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        for (int site = 0; site < GetNsite(); site++) {
                 delete[] PCAT[site];
        }

        delete[] PCAT;
	delete[] tmpVI;
	delete[] tmp;
}


double PoissonSBDPProfileProcessVI::MoveProCluster() {
        
        // UpdateOccupancyNumbers();
        ResampleWeights();
		int site = (int) (GetNsite() * rnd::GetRandom().Uniform());

		int bk = alloc[site];
		int h = occupancy[bk] > 1 ? Ncomponent+1 : Ncomponent;

		// make a new mode Nmode <= i < h
		/*for (int i=Ncomponent; i<h ; i++)	{
			CreateComponent(i);
		}

		RemoveSite(site,bk);*/

		// Gibbs

		double total = 0;
		double* mModeGibbsGrid = new double[h];
		double* mLogSamplingArray = new double[GetNmodeMax()];

		double max = 0;
        for (int site = 0; site < GetNsite(); site++) {
		/*for (int mode = 0; mode < GetNmodeMax(); mode++)	{
			mLogSamplingArray[mode] = LogStatProb(site,mode);
                        if ((!mode) || (max < mLogSamplingArray[mode])) {
                                      max = mLogSamplingArray[mode];
                        }
		}*/
		for (int mode = 0; mode < GetNmodeMax(); mode++)	{
			PCAT[site][mode] = exp(weightVI[mode] + LogStatProb(site,mode));
                        // PCAT[site][mode] = rnd::GetRandom().Uniform();
		}
                /*double mean = 0;
                for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                       mean += PCAT[site][mode];
                }
                mean /= GetNmodeMax();

                maxPCAT = 0;
                for (int mode = 0; mode < GetNmodeMax(); mode++)	{
                        if ((!mode) || (maxPCAT < PCAT[site][mode]))  {
                               maxPCAT = PCAT[site][mode];
                        }
                }
                double minPCAT = 1;
                for (int mode = 0; mode < GetNmodeMax(); mode++)	{
                        if ((!mode) || (minPCAT > PCAT[site][mode]))  {
                               minPCAT = PCAT[site][mode];
                        }
                }*/

                /*double totalsquare = 0;
                       for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                             totalsquare += (PCAT[site][mode] * PCAT[site][mode]);
                       }
                double norm = sqrt(totalsquare);*/

                double total = 0;
                       for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                             total += PCAT[site][mode];
                       }
                double norm = total;

                for (int mode = 0; mode < GetNmodeMax(); mode++)	{
                        NormalizePCAT[site][mode] = PCAT[site][mode] / norm;
                        // NormalizePCAT[site][mode] = (PCAT[site][mode] - minPCAT) / (maxPCAT - minPCAT);
                        // NormalizePCAT[site][mode] = (PCAT[site][mode]) / (PCAT[site][mode]-1);
                }
                double MaxProCluster = 0;
                for (int mode = 0; mode < GetNmodeMax(); mode++)	{
                        if ((!mode) || (MaxProCluster < NormalizePCAT[site][mode]))  {
                               MaxProCluster = NormalizePCAT[site][mode];
                               alloc[site] = mode;
                        }
                }
        }
		/*if (mode >= Ncomponent)	{
			if (mode > Ncomponent)	{
				SwapComponents(mode, Ncomponent);
			}
			Ncomponent++;
		}
		if (! occupancy[bk])	{
			if (bk != Ncomponent-1)	{
				SwapComponents(bk,Ncomponent-1);
			}
			Ncomponent--;
		}*/

		delete[] mModeGibbsGrid;
		delete[] mLogSamplingArray;

	return 1.0;

}




































