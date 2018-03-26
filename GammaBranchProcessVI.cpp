
/********************



**********************/


#include "GammaBranchProcessVI.h"
#include <iostream>
#include <cassert>
#include "Parallel.h"
#include <string.h>

void GammaBranchProcessVI::ToStream(ostream& os)	{

	SetNamesFromLengths();
	tree->ToStream(os);
	os << branchalpha << '\n';
	os << branchbeta << '\n';
	/*
	for (int j=0; j<GetNbranch(); j++)	{
		os << blarray[j] << '\t';
	}
	os << '\n';
	*/
}

void GammaBranchProcessVI::ToStreamWithLengths(ostream& os, const Link* from)	{

	if (from->isLeaf())	{
		os << from->GetNode()->GetName();
	}
	else	{
		os << "(";
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStreamWithLengths(os, link->Out());
			if (link->Next() != from)	{
				os << ",";
			}
		}
		os << ")";
	}
	if (! from->isRoot())	{
		os << ":" << blarray[from->GetBranch()->GetIndex()];
	}
}


void GammaBranchProcessVI::FromStream(istream& is)	{

	tree->ReadFromStream(is);
	tree->RegisterWith(tree->GetTaxonSet());
	SetLengthsFromNames();
	branchalpha = -1;
	branchbeta = -1;
	is >> branchalpha;
	is >> branchbeta;
	/*
	for (int j=0; j<GetNbranch(); j++)	{
		is >> blarray[j];
	}
	*/
}

// -------------- Compute the variational distribution of GammaBranchProcess ------------------------------------------------------------------------------------ 
	
double GammaBranchProcessVI::LogBranchLengthPrior(const Branch* branch)	{
	int index = branch->GetIndex();
	return branchalphaVI[index] * log(branchbetaVI[index]) - rnd::GetRandom().logGamma(branchalphaVI[index]) + (branchalphaVI[index]-1) * log(blarray[index]) - branchbetaVI[index] * blarray[index];
}

double GammaBranchProcessVI::ELBOLength()  {
 
        ELBO_Length = 0;  
        for (int i=1; i<GetNbranch(); i++)   {
             ELBO_Length += GetBranchLengthSuffStatCount(i) * ((gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))));
             ELBO_Length += (branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha - 1) * (gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))) - branchbeta * (GetbranchalphaVI(i) / GetbranchbetaVI(i)));
             ELBO_Length -= (GetbranchalphaVI(i) * log(GetbranchbetaVI(i)) - rnd::GetRandom().logGamma(GetbranchalphaVI(i)));
             ELBO_Length -= ((GetbranchalphaVI(i) - 1) * (gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))) - GetbranchbetaVI(i) * (GetbranchalphaVI(i) / GetbranchbetaVI(i)));
        }
        return ELBO_Length;
}

/*void GammaBranchProcessVI::SampleLength(const Branch* branch)	{
	int index = branch->GetIndex();
	blarray[index] = rnd::GetRandom().Gamma(branchalpha,branchbeta);
}*/

void GammaBranchProcessVI::SampleLength(const Branch* branch)   {

        int index = branch->GetIndex();
        blarray[index] = rnd::GetRandom().Gamma(branchalphaVI[index], branchbetaVI[index]);
}
	
void GammaBranchProcessVI::SampleLength()	{
	cerr << "sample length\n";
	exit(1);
	branchalpha = rnd::GetRandom().sExpo();
	branchbeta = rnd::GetRandom().sExpo();
	blarray[0] = 0;
	// SampleLength();
}

/*void GammaBranchProcessVI::SampleLength()  {
       branchalpha = 1;
       branchbeta = 10;
       for (int i=1; i<GetNbranch(); i++)  {
		blarray[i] = rnd::GetRandom().Gamma(branchalpha, branchbeta);
       }
}*/


// ------------- Compute the derivative of variational distribution ------------------------------------------------------------------------------------------------------- 

double GammaBranchProcessVI::DQBranchAlpha(const Branch* branch) {
         int index = branch->GetIndex();
         return dqbranchalpha = -gsl_sf_psi(branchalpha) + log(branchbeta) + log(blarray[index]);
}
 
double GammaBranchProcessVI::DQBranchBeta(const Branch* branch) {
         int index = branch->GetIndex();
         return branchalpha/branchbeta - blarray[index];
}

//---------------------------------------------------------------------------------------------------------------------------
//                             Stochastic Optimization for Variational Parameters  
//---------------------------------------------------------------------------------------------------------------------------

void GammaBranchProcessVI::TotalParameterLength()  {

      TotalParaLength = 0;
      for(int i=1; i<GetNbranch(); i++) {
           TotalParaLength += (branchalphaVI[i] / branchbetaVI[i]);
      }
}

double GammaBranchProcessVI::MoveLength()	{

	for (int i=1; i<GetNbranch(); i++)	{
		blarray[i] = rnd::GetRandom().Gamma(GetbranchalphaVI(i), GetbranchbetaVI(i));
	}
        return 1.0;
}

void GammaBranchProcessVI::EstimateBranchParameter()  {
     
     // GlobalUpdateBranchLengthSuffStat();
     // TotalParameterRate();
     for (int i=1; i<GetNbranch(); i++)	{ 
     branchalphahat[i] = branchalpha + GetBranchLengthSuffStatCount(i);
     branchbetahat[i] = branchbeta + GetBranchLengthSuffStatBeta(i);
     }  
}

double GammaBranchProcessVI::MoveBranchParamter()  {
       for (int i=1; i<GetNbranch(); i++) { 
       
       LearningRateBranch[i] = pow(taubranch + i, -1 * kappabranch); 
       branchalphaVI[i] = ((1 - LearningRateBranch[i]) * branchalphaVI[i]) + (LearningRateBranch[i] * Getbranchalphahat(i));
       branchbetaVI[i] = ((1 - LearningRateBranch[i]) * branchbetaVI[i]) + (LearningRateBranch[i] * Getbranchbetahat(i)); 
       
       /*branchalphaVI[i] = 1;
       branchbetaVI[i] = 10;*/
       }
   return 1.0;
}

// ---------------------------------------------------------------------------------------------------------------------------
//             Parallel Computation (Open-MPI) for Stochastic Optimization for Variational Parameters  
// ---------------------------------------------------------------------------------------------------------------------------

double GammaBranchProcessVI::Move()  {

       // MoveLength();
       // EstimateBranchParameter();  
       // MoveBranchParamter();
       // MoveLength();
       // GlobalMoveLength();
       // ELBOLength();  
       GlobalMoveLength();
       GlobalUpdateBranchLengthSuffStat();
       GlobalEstimateBranchParameter();
       GlobalMoveBranchParamter();
       GlobalELBOLength();
       // ELBOLength();
       return 1.0;
}

void GammaBranchProcessVI::SlaveELBOLength() {

       assert(GetMyid() > 0);
       int brwidth = GetNbranch()/(GetNprocs()-1);
       int brmin[GetNprocs()-1];
       int brmax[GetNprocs()-1];
       for (int i=0; i<GetNprocs()-1; ++i)  {
               brmin[i] = brwidth*i;
               brmax[i] = brwidth*(1+i);
               if (i == (GetNprocs()-2)) brmax[i] = GetNbranch();
       }  

       double totalelbo = 0;
       for (int i=brmin[GetMyid()-1]; i<brmax[GetMyid()-1]; i++)   {
             totalelbo += GetBranchLengthSuffStatCount(i) * ((gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))));
             totalelbo += (branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha - 1) * (gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))) - branchbeta * (GetbranchalphaVI(i) / GetbranchbetaVI(i)));
             totalelbo -= (GetbranchalphaVI(i) * log(GetbranchbetaVI(i)) - rnd::GetRandom().logGamma(GetbranchalphaVI(i)));
             totalelbo -= ((GetbranchalphaVI(i) - 1) * (gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))) - GetbranchbetaVI(i) * (GetbranchalphaVI(i) / GetbranchbetaVI(i)));
        }
        MPI_Send(&totalelbo,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

/*void GammaBranchProcessVI::SlaveELBOLength()     {

       assert(GetMyid() > 0);
       int brwidth = GetNbranch()/(GetNprocs()-1);
       int brmin[GetNprocs()-1];
       int brmax[GetNprocs()-1];
       for (int i=0; i<GetNprocs()-1; i++)     {
              brmin[i] = brwidth * i;
              brmax[i] = brwidth * (1+i);
              if (i == (GetNprocs()-2)) brmax[i] = GetNbranch();
       }
       double totsuffstat = 0;
       for (int i=0; i<GetNbranch(); i++)    {
             totsuffstat += GetBranchLengthSuffStatCount(i);
       }
       double totalelbo = 0;
       for (int i=brmin[GetMyid()-1]; i<brmax[GetMyid()-1]; i++)     {
             double branchsuff = GetBranchLengthSuffStatCount(i) / totsuffstat;
             totalelbo += branchsuff * ((gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))));
             totalelbo += (branchalpha * log(branchbeta) - rnd::GetRandom().logGamma(branchalpha) + (branchalpha - 1) * (gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))) - branchbeta * (GetbranchalphaVI(i) / GetbranchbetaVI(i)));
             totalelbo -= (GetbranchalphaVI(i) * log(GetbranchbetaVI(i)) - rnd::GetRandom().logGamma(GetbranchalphaVI(i)));
             totalelbo -= ((GetbranchalphaVI(i) - 1) * (gsl_sf_psi(GetbranchalphaVI(i)) - log(GetbranchbetaVI(i))) - GetbranchbetaVI(i) * (GetbranchalphaVI(i) / GetbranchbetaVI(i)));
       }
       MPI_Send(&totalelbo,1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}*/

double GammaBranchProcessVI::GlobalELBOLength() {

       assert(GetMyid() == 0);
       MPI_Status stat;
       MESSAGE signal = ELBOLengthVI;
       
       MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
       int brwidth = GetNbranch()/(GetNprocs()-1);
       int brmin[GetNprocs()-1];
       int brmax[GetNprocs()-1];
       for (int i=0; i<GetNprocs()-1; ++i)  {
                 brmin[i] = brwidth*i;
                 brmax[i] = brwidth*(1+i);
                 if (i == (GetNprocs() - 2)) brmax[i] = GetNbranch();
       }

       ELBO_Length = 0.0;
       double sum;
       for (int i=1; i<GetNprocs(); ++i)   {
               MPI_Recv(&sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
               ELBO_Length += sum;
       }
       return ELBO_Length;
}

void GammaBranchProcessVI::GlobalMoveLength()  {
       assert(GetMyid() == 0);
       MPI_Status stat;
       MESSAGE signal = UPDATE_MOVEBRAN;

       MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
       // split Nbranch() among GetNprocs()-1 slaves
       int brwidth = GetNbranch()/(GetNprocs()-1);
       int brmin[GetNprocs()-1];
       int brmax[GetNprocs()-1];
       for(int i=0; i<GetNprocs()-1; ++i)  {
                brmin[i] = brwidth*i;
                brmax[i] = brwidth*(1+i);
                if (i == (GetNprocs()-2)) brmax[i] = GetNbranch();
       }
       double* tmpbran = new double[GetNbranch() + 1];
       for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpbran,(brmax[i-1] - brmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int k = 0;
               for(int j=brmin[i-1]; j<brmax[i-1]; ++j) {
                       blarray[j] = tmpbran[k];
                       k++;
               }
       }
       delete[] tmpbran;
}

void GammaBranchProcessVI::SlaveMoveLength()  {
        assert(GetMyid() > 0);
        double* tmpbran = new double[GetNbranch() + 1];
        // split Nbranch among GetNprocs() -1 slaves 
       int brwidth = GetNbranch()/(GetNprocs()-1);
       int brmin[GetNprocs()-1];
       int brmax[GetNprocs()-1];
       for(int i=0; i<GetNprocs()-1; ++i) {
              brmin[i] = brwidth*i;
              brmax[i] = brwidth*(1+i);
              if (i == (GetNprocs()-2)) brmax[i] = GetNbranch();
       }
       int k = 0;
       for(int j=brmin[GetMyid()-1]; j<brmax[GetMyid()-1]; j++) {
              blarray[j] = rnd::GetRandom().Gamma(GetbranchalphaVI(j), GetbranchbetaVI(j));
              tmpbran[k] = blarray[j];
              k++;
       }
       MPI_Send(tmpbran,(brmax[GetMyid()-1] - brmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       delete[] tmpbran;
}

void GammaBranchProcessVI::GlobalEstimateBranchParameter()	{

	assert(GetMyid() == 0);
	MPI_Status stat;
	MESSAGE signal = UPDATE_ESTBRANCH;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
        // split Nbranch() among GetNprocs()-1 slaves
        int brwidth = GetNbranch()/(GetNprocs()-1);
        int brmin[GetNprocs()-1];
        int brmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i)  {
                brmin[i] = brwidth*i;
                brmax[i] = brwidth*(1+i);
                if (i == (GetNprocs()-2)) brmax[i] = GetNbranch();
        }
        double* tmpbal = new double[GetNbranch() + 1];
        double* tmpbbe = new double[GetNbranch() + 1];
        for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpbal,(brmax[i-1] - brmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int k = 0;
               for(int j=brmin[i-1]; j<brmax[i-1]; ++j) {
                       branchalphahat[j] = tmpbal[k];
                       k++;
               }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpbbe,(brmax[i-1] - brmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int l = 0;
               for(int j=brmin[i-1]; j<brmax[i-1]; ++j) {
                       branchbetahat[j] = tmpbbe[l];
                       l++;
               }
        }
        delete[] tmpbal;
        delete[] tmpbbe;
}

void GammaBranchProcessVI::SlaveEstimateBranchParameter()	{
     
       assert(GetMyid() > 0);
       double* tmpbal = new double[GetNbranch() + 1];
       double* tmpbbe = new double[GetNbranch() + 1];

       // split Nbranch among GetNprocs() -1 slaves 
       int brwidth = GetNbranch()/(GetNprocs()-1);
       int brmin[GetNprocs()-1];
       int brmax[GetNprocs()-1];
       for(int i=0; i<GetNprocs()-1; ++i) {
              brmin[i] = brwidth*i;
              brmax[i] = brwidth*(1+i);
              if (i == (GetNprocs()-2)) brmax[i] = GetNbranch();
       }
       int k = 0;
       for(int j=brmin[GetMyid()-1]; j<brmax[GetMyid()-1]; j++) {
              branchalphahat[j] = branchalpha + GetBranchLengthSuffStatCount(j);
              tmpbal[k] = branchalphahat[j];
              k++;
       }
       int l = 0;
       for(int j=brmin[GetMyid()-1]; j<brmax[GetMyid()-1]; j++) {
              branchbetahat[j] = branchbeta + GetBranchLengthSuffStatBeta(j);
              tmpbbe[l] = branchbetahat[j];
              l++; 
       }

       MPI_Send(tmpbal,(brmax[GetMyid()-1] - brmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       MPI_Barrier(MPI_COMM_WORLD);
       MPI_Send(tmpbbe,(brmax[GetMyid()-1] - brmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       delete[] tmpbal;
       delete[] tmpbbe;
}

void GammaBranchProcessVI::GlobalMoveBranchParamter()  {

       assert(GetMyid() == 0);
       MPI_Status stat;
       MESSAGE signal = UPDATE_MOVEBRANCH;

       MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
       // split Nbranch() among GetNprocs()-1 slaves
       int movwidth = GetNbranch()/(GetNprocs()-1);
       int movmin[GetNprocs()-1];
       int movmax[GetNprocs()-1];
        for(int i=0; i<GetNprocs()-1; ++i)  {
                movmin[i] = movwidth*i;
                movmax[i] = movwidth*(1+i);
                if (i == (GetNprocs()-2)) movmax[i] = GetNbranch();
        }
        double* tmpmoval = new double[GetNbranch() + 1];
        double* tmpmovbe = new double[GetNbranch() + 1];
        for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpmoval,(movmax[i-1] - movmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int k = 0;
               for(int j=movmin[i-1]; j<movmax[i-1]; ++j) {
                       branchalphaVI[j] = tmpmoval[k];
                       k++;
               }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=1; i<GetNprocs(); ++i)  {
               MPI_Recv(tmpmovbe,(movmax[i-1] - movmin[i-1])+1,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
               int l = 0;
               for(int j=movmin[i-1]; j<movmax[i-1]; ++j) {
                       branchbetaVI[j] = tmpmovbe[l];
                       l++;
               }
        }
        delete[] tmpmoval;
        delete[] tmpmovbe;
}

void GammaBranchProcessVI::SlaveMoveBranchParamter()  {

       assert(GetMyid() > 0);
       double* tmpmoval = new double[GetNbranch() + 1];
       double* tmpmovbe = new double[GetNbranch() + 1];

       // Split Nbranch among GetNprocs() - 1 slaves
       int movwidth = GetNbranch()/(GetNprocs()-1);
       int movmin[GetNprocs()-1];
       int movmax[GetNprocs()-1];
       for(int i=0; i<GetNprocs()-1; ++i) {
             movmin[i] = movwidth*i;
             movmax[i] = movwidth*(1+i);
             if (i == (GetNprocs()-2)) movmax[i] = GetNbranch();
       }
       for(int j=movmin[GetMyid()-1]; j<movmax[GetMyid()-1]; j++) {
             LearningRateBranch[j] = pow(taubranch + j, -1 * kappabranch); 
       }
       int k = 0;
       for(int j=movmin[GetMyid()-1]; j<movmax[GetMyid()-1]; j++) { 
              branchalphaVI[j] = ((1 - LearningRateBranch[j]) * branchalphaVI[j]) + (LearningRateBranch[j] * Getbranchalphahat(j)); 
              tmpmoval[k] = branchalphaVI[j];
              k++;
       }
       int l = 0;
       for(int j=movmin[GetMyid()-1]; j<movmax[GetMyid()-1]; j++) {
              branchbetaVI[j] = ((1 - LearningRateBranch[j]) * branchbetaVI[j]) + (LearningRateBranch[j] * Getbranchbetahat(j)); 
              tmpmovbe[l] = branchbetaVI[j];
              l++;
       } 

       MPI_Send(tmpmoval,(movmax[GetMyid()-1] - movmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       MPI_Barrier(MPI_COMM_WORLD);
       MPI_Send(tmpmovbe,(movmax[GetMyid()-1] - movmin[GetMyid()-1])+1,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
       delete[] tmpmoval;
       delete[] tmpmovbe;
}


/*void GammaBranchProcessVI::GlobalMoveLength()  {

       assert(GetMyid() == 0);
       int i,j,nbranch = GetNbranch(),nprocs = GetNprocs();
       MPI_Status stat;
       MESSAGE signal = UPDATE_MOVEBRAN;

       MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);
       
       #ifdef BYPE_COM
       double k,l;
       unsigned char* bvector = new unsigned char[nbranch*(sizeof(double))];
       for(i=1; i<nprocs; ++i) {
               MPI_Recv(bvector,nbranch*(sizeof(double)),MPI_UNSIGNED_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
               for(j=0; i<nbranch; ++j) {
                     l = 0;
                     for(k=sizeof(double)-1; k>=0; --k) {
                             l = (l << 8) + bvector[sizeof(double)*j+k];
                     }
                     blarray[j] = l;
               }
       }
       delete[] bvector;
       #else
       double ivector[nbranch];
       for(i =1; i<nprocs; ++i) {
                MPI_Recv(ivector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
                for(j=0; j<nbranch; ++j) {
                        blarray[j] = ivector[j]; 
                }
       }
       #endif
       
}

void GammaBranchProcessVI::SlaveMoveLength()  {

        assert(GetMyid() > 0);

        MoveLength();

        int workload = GetNbranch();
        #ifdef BYTE_COM
        int i,n = 0;
        unsigned int j;
        unsigned char el_int[sizeof(double)];
        unsigned char* bvector = new unsigned char[workload*(sizeof(double))];
        for(i=0; i<workload; ++i) {
                convert(el_int,blarray[i]);
                for(j =0; j<sizeof(double); ++j) {
                       bvector[n] = el_int[j]; n++;
                } 
        }
        MPI_Send(bvector,workload*(sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
        delete[] bvector;
        #else
        MPI_Send(blarray,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
        #endif
}

void GammaBranchProcessVI::GlobalMoveBranchParamter()	{

	assert(GetMyid() == 0);
	int i,j,nbranch = GetNbranch(),nprocs = GetNprocs();
	MPI_Status stat;
	MESSAGE signal = UPDATE_MOVEBRANCH;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	#ifdef BYTE_COM
	// should be summed over all slaves (reduced)
	double k,l;
	double x;
	unsigned char* bvector = new unsigned char[nbranch*(sizeof(double)+sizeof(double))];

	for(i=1; i<nprocs; ++i) {
		MPI_Recv(bvector,nbranch*(sizeof(double)+sizeof(double)),MPI_UNSIGNED_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			l = 0;
			for(k=sizeof(double)-1; k>=0; --k) {
				l = (l << 8) + bvector[sizeof(double)*j+k]; 
			}
			branchalphaVI[j] = l;
		}
		for(j=0; j<nbranch; ++j) {
			memcpy(&x,&bvector[sizeof(double)*nbranch+sizeof(double)*j],sizeof(double));
			branchbetaVI[j] = x;
		}
	}
	delete[] bvector;
	#else
	double ivector[nbranch];
	double dvector[nbranch];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(ivector,nbranch,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchalphaVI[j] = ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(dvector,nbranch,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchbetaVI[j] = dvector[j];
		}
	}
	
	#endif
}

void GammaBranchProcessVI::SlaveMoveBranchParamter()	{
        
        assert(GetMyid() > 0);

	MoveBranchParamter();
	
	int workload = GetNbranch();
	#ifdef BYTE_COM
	int i,n = 0;
	unsigned int j;
	unsigned char el_int[sizeof(double)],el_dbl[sizeof(double)];
	unsigned char* bvector = new unsigned char[workload*(sizeof(double)+sizeof(double))];
	for(i=0; i<workload; ++i) {
		convert(el_int,branchalphaVI[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_int[j]; n++;
		}
	}
	for(i=0; i<workload; ++i) {
		convert(el_dbl,branchbetaVI[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_dbl[j]; n++;
		}
	}
	MPI_Send(bvector,workload*(sizeof(double)+sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;
	#else
	MPI_Send(branchalphaVI,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(branchbetaVI,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	#endif
}

void GammaBranchProcessVI::GlobalEstimateBranchParameter()	{

	assert(GetMyid() == 0);
	int i,j,nbranch = GetNbranch(),nprocs = GetNprocs();
	MPI_Status stat;
	MESSAGE signal = UPDATE_ESTBRANCH;

	MPI_Bcast(&signal,1,MPI_INT,0,MPI_COMM_WORLD);

	#ifdef BYTE_COM
	// should be summed over all slaves (reduced)
	double k;
        double l;
	double x;
	unsigned char* bvector = new unsigned char[nbranch*(sizeof(double)+sizeof(double))];

	for(i=1; i<nprocs; ++i) {
		MPI_Recv(bvector,nbranch*(sizeof(double)+sizeof(double)),MPI_UNSIGNED_CHAR,i,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			l = 0;
			for(k=sizeof(double)-1; k>=0; --k) {
				l = (l << 8) + bvector[sizeof(double)*j+k]; 
			}
			branchalphahat[j] = l;
		}
		for(j=0; j<nbranch; ++j) {
			memcpy(&x,&bvector[sizeof(double)*nbranch+sizeof(double)*j],sizeof(double));
			branchbetahat[j] = x;
		}
	}
	delete[] bvector;
	#else
	double ivector[nbranch];
	double dvector[nbranch];
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(ivector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(ivector,nbranch,MPI_INT,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchalphahat[j] = ivector[j];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(i=1; i<nprocs; ++i) {
		MPI_Recv(dvector,nbranch,MPI_DOUBLE,i,TAG1,MPI_COMM_WORLD,&stat);
		// MPI_Recv(dvector,nbranch,MPI_DOUBLE,MPI_ANY_SOURCE,TAG1,MPI_COMM_WORLD,&stat);
		for(j=0; j<nbranch; ++j) {
			branchbetahat[j] = dvector[j];
		}
	}
	
	#endif
}

void GammaBranchProcessVI::SlaveEstimateBranchParameter()	{
        
        assert(GetMyid() > 0);

	EstimateBranchParameter();
	
	int workload = GetNbranch();
	#ifdef BYTE_COM
	int i,n = 0;
	unsigned int j;
	unsigned char el_int[sizeof(double)],el_dbl[sizeof(double)];
	unsigned char* bvector = new unsigned char[workload*(sizeof(double)+sizeof(double))];
	for(i=0; i<workload; ++i) {
		convert(el_int,branchalphahat[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_int[j]; n++;
		}
	}
	for(i=0; i<workload; ++i) {
		convert(el_dbl,branchbetahat[i]);
		for(j=0; j<sizeof(double); ++j) {
			bvector[n] = el_dbl[j]; n++;
		}
	}
	MPI_Send(bvector,workload*(sizeof(double)+sizeof(double)),MPI_UNSIGNED_CHAR,0,TAG1,MPI_COMM_WORLD);
	delete[] bvector;
	#else
	MPI_Send(branchalphahat,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Send(branchbetahat,workload,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);

	#endif
}*/

// -------------------------------------------------------------------------------------------------------------------------------------------------------



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/* double GammaBranchProcess::Move(double tuning, int nrep)	{
	double total = MoveLength();
	total += MoveBranchBeta(tuning,nrep);
	// cerr << "in bl move : " << GetTotalLength() << '\n';
	return total;
}

double GammaBranchProcess::MoveBranchBeta(double tuning, int nrep)	{
	int Naccepted = 0;
	for (int rep=0; rep<nrep; rep++)	{
		double deltalogprob = - LogHyperPrior() - LogLengthPrior() - LengthSuffStatLogProb();
		double m = tuning * (rnd::GetRandom().Uniform() - 0.5);
		double e = exp(m);
		branchbeta *= e;
		deltalogprob += LogHyperPrior() + LogLengthPrior() + LengthSuffStatLogProb();
		deltalogprob += m;
		int accepted = (log(rnd::GetRandom().Uniform()) < deltalogprob);
		if (accepted)	{
			Naccepted ++;
		}
		else	{
			branchbeta /= e;
		}
	}
	return ((double) Naccepted) / nrep;
}



*/
