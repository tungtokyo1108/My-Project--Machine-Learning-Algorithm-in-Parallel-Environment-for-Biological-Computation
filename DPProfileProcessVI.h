
/********************

Copyright 2016-2018 Dang Thanh Tung, Hirohisa Kishino
Fast Computational version of PhyloBayes MPI is free software which replace the sampling algorithm by fast optimization  
algorithm : you can redistribute it and/or modify it under the terms of the GNU General Public License.

**********************/


#ifndef DPPROFILEVI_H
#define DPPROFILEVI_H

#include <cmath>
#include "MixtureProfileProcessVI.h"

// general superclass for all finite process mixtures on site-specific profiles
class DPProfileProcessVI: public virtual MixtureProfileProcessVI	{

	public:

	DPProfileProcessVI() : kappa(1), movekappa(true), kappaprior(0) {}
	virtual ~DPProfileProcessVI(){}

	protected:

	virtual double IncrementalDPMove(int nrep) = 0;
	double MoveHyper(double tuning, int nrep);
	virtual double MoveKappa(double tuning, int nrep);


	// static allocation of many component-specific variables
	// such as: profiles, occupancy number
	// basically everything except substitution matrices

	// called at the beginning and end of the run (see PhyloProcess)
	/*
	virtual void Create(int innsite, int indim);
	virtual void Delete();
	*/

	// multinomial 
	virtual double LogProxy(int site, int cat) = 0;
	virtual void SampleAlloc();
	void SampleHyper();

	// kappa has an exponential prior of mean 10
	double LogHyperPrior();
	virtual double LogAllocPrior();

	virtual void DrawProfileFromPrior()	{
		cerr << "error: in DPProfileProcess::DrawProfileFromPrior\n";
		exit(1);
	}

	double kappa;
        /*double* kappa_alpha;
        double* kappa_beta;*/

	bool movekappa;
	int kappaprior;
	// 0 : exponential of mean 20
	// 1 : jeffreys prior 
	int dirweightprior;
	// 0 : flexible
	// 1 : rigid
};

#endif

