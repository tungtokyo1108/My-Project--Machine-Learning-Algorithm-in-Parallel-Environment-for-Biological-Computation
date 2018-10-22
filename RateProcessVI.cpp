
/********************

**********************/


#include "RateProcessVI.h"
#include <cassert>
#include "Parallel.h"
#include <string.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* RateProcess
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

double RateProcessVI::GetMeanRate()	{
	double total = 0;
	for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{
	// for (int i=0; i<GetNsite(); i++)	{
		total += GetRate(i);
	}
	return total / (GetSiteMax() - GetSiteMin());
	// return total / GetNsite();
}

double RateProcessVI::GetMeanRateVI()	{
       double total = 0;
       // for(int i=0; i<GetNsite(); i++) {
       for (int i=GetSiteMin(); i<GetSiteMax(); i++)	{ 
               total += (GetRateAlphaVI(i)/GetRateBetaVI(i));
       }
       return total;
}
