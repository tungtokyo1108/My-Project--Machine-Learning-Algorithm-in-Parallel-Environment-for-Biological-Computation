
/********************

Adopted PhyloBayes MPI. https://github.com/bayesiancook/pbmpi
Lartillot, N., Rodrigue, N., Stubbs, D. & Richer, J. 
PhyloBayes MPI: Phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Syst. Biol. 62, 611–615 (2013).

**********************/

#ifndef CHRONO_H
#define CHRONO_H

#include <ctime>
#include <sys/time.h>
#include <iostream>
#include <cstdlib>

using namespace std;

class Chrono	{

	public:

	Chrono(){
		Reset();
	};
	~Chrono() {};
	void Reset();
	void Start();
	void Stop();

	void ToStream(ostream& os);
	void FromStream(istream& is);

	int operator++();

	double GetTime();
	double GetTimePerCount();
	int GetCount();

	private:

	// this is in milli seconds
	double sec1;
	double sec2;
	double milli1;
	double milli2;

	double TotalTime;
	int N;

}
;


#endif		
