
/********************


**********************/

#ifndef POISSONDPPROFILEVI_H
#define POISSONDPPROFILEVI_H

#include "PoissonMixtureProfileProcessVI.h"
#include "DPProfileProcessVI.h"

// superclass for Poisson (F81) implementations
class PoissonDPProfileProcessVI: public virtual PoissonMixtureProfileProcessVI, public virtual DPProfileProcessVI	{

	public:

	PoissonDPProfileProcessVI() {}
	virtual ~PoissonDPProfileProcessVI() {}

	virtual double Move(double tuning = 1, int n = 1, int nrep = 1)	{
		for (int rep=0; rep<nrep; rep++)	{
			GlobalUpdateParameters();
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			//IncrementalDPMove(2);
			//MoveProfile();
			//MoveHyper(tuning,10);
			//MoveHyper(0.1 * tuning,10);
		}
		return 1;
	/*
			GlobalUpdateSiteProfileSuffStat();
			// UpdateSiteProfileSuffStat();
			// GlobalUpdateParameters();
			// GlobalUpdateSiteProfileSuffStat();
			UpdateModeProfileSuffStat();
			IncrementalDPMove(nrep);
			MoveProfile();
			MoveHyper(tuning,nrep);
	*/
	}

	virtual void ToStream(ostream& os);
	virtual void FromStream(istream& is);

	protected:

	
	virtual void Create(int innsite, int indim)	{
		DPProfileProcessVI::Create(innsite,indim);
		PoissonMixtureProfileProcessVI::Create(innsite,indim);
	}

	virtual void Delete()	{
		DPProfileProcessVI::Delete();
		PoissonMixtureProfileProcessVI::Delete();
	}
	

	//double IncrementalDPMove(int nrep);

};

#endif

