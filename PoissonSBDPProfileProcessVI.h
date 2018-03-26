
/********************


**********************/

#ifndef POISSONSBDPPROFILEVI_H
#define POISSONSBDPPROFILEVI_H

#include "PoissonDPProfileProcessVI.h"
#include "SBDPProfileProcessVI.h"
#include "PoissonMixtureProfileProcessVI.h"
#include "Random.h"
#include <gsl/gsl_sf_psi.h>

// superclass for Poisson (F81) implementations
class PoissonSBDPProfileProcessVI: public virtual PoissonDPProfileProcessVI, public virtual SBDPProfileProcessVI,  public virtual PoissonMixtureProfileProcessVI {

	public:

	PoissonSBDPProfileProcessVI() : InitIncremental(0) {}
	virtual ~PoissonSBDPProfileProcessVI() {}
        double maxPCAT;
        double GetmaxPCAT() {return maxPCAT;}
        double TotalPCAT() {
             double total = 0;
                for (int site = 0; site < 1; site++) {
                      for (int mode = 0; mode < GetNmodeMax(); mode++)    {
                           // total += (NormalizePCAT[site][mode] * NormalizePCAT[site][mode]);
                           total += NormalizePCAT[site][mode];
                      }
                }
            return sqrt(total);
        }
	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{

		// totchrono.Start();
		for (int rep=0; rep<nrep; rep++)	{

			// incchrono.Start();
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			if ((!rep) && InitIncremental)	{
				cerr << "init incremental\n";
				InitIncremental--;
				IncrementalSampleAlloc();
				UpdateModeProfileSuffStat();
			}
                        // MoveProCluster(); 
			GlobalMixMove(1,1,0.001);
			// MoveOccupiedCompAlloc(5);
			// MoveAdjacentCompAlloc(5);
			// incchrono.Stop();
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// MoveHyper(tuning,1);
                        SBDPProfileProcessVI::MoveKappa(1,1);
                        PoissonMixtureProfileProcessVI::MoveDirWeightVI(); 
		}
		// totchrono.Stop();
		return 1;
	}

	virtual double LogProxy(int site, int cat)	{
		return DiffLogSampling(cat,site);
	}
        double** PCAT;
        // double** TransPCAT;
        double** NormalizePCAT;
        double* UpdateProCAT;
        double GlobalMixMove(int nrep, int nallocrep, double epsilon);
	void SlaveMixMove();
        double GetPCAT(int cat) {return UpdateProCAT[cat];}
        double Getproclus()   {
               double total = 0;
               for (int mode = 0; mode < GetNmodeMax(); mode++)    {
                    total += UpdateProCAT[mode];
               }
               return total;
        }
        double GetTotalSitePCAT()  {

                double** TransPCAT = new double*[GetNmodeMax()];
                for (int j = 0; j < GetNmodeMax(); j++)  {
                         TransPCAT[j] = new double[GetNsite()];
                }
 
                for (int site = 0; site < GetNsite(); site++) {
                      for (int mode = 0; mode < GetNmodeMax(); mode++)    {
                            TransPCAT[mode][site] = NormalizePCAT[site][mode];
                      }
                }
              
                for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                           UpdateProCAT[mode] = 0;
                } 
                for (int mode = 0; mode < GetNmodeMax(); mode++)  {
                       for (int site = 0; site < GetNsite(); site++)  {
                            UpdateProCAT[mode] += TransPCAT[mode][site];
                       }
                } 

                for (int j=0; j<GetNmodeMax(); j++) {
                         delete[] TransPCAT[j];
                }
                delete[] TransPCAT;

            return 1.0;     
        }
 
        double MoveProCluster();
	protected:

	virtual void Create(int innsite, int indim);
		/*PoissonDPProfileProcessVI::Create(innsite,indim);
		SBDPProfileProcessVI::Create(innsite,indim);
                
                PCAT = new double[GetNmodeMax()];
                for (int mode = 0; mode < GetNmodeMax(); mode++) {
                     
                     PCAT[mode] = rnd::GetRandom().Uniform();
                     // PCAT[mode] = 0.5;
                     // PCAT[mode] = exp(GetProfileParameter());
                }*/ 
	

	virtual void Delete()	{
		SBDPProfileProcessVI::Delete();
		PoissonDPProfileProcessVI::Delete();
	}

	/*double GlobalMixMove(int nrep, int nallocrep, double epsilon);
	void SlaveMixMove();*/

	double IncrementalDPMove(int nrep, double c)	{
		cerr << "error : in poisson sbdp incremental\n";
		exit(1);
		return 1;
	}
	double IncrementalDPMove(int nrep)	{
		cerr << "error : in poisson sbdp incremental\n";
		exit(1);
		return 1;
	}

	virtual void SwapComponents(int cat1, int cat2)	{
		SBDPProfileProcessVI::SwapComponents(cat1,cat2);
	}

	// virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is)	{
		PoissonDPProfileProcessVI::FromStream(is);
		ResampleWeights();
	}

	int InitIncremental;
        
        
        
        

};

#endif

