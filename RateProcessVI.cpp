
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

/*double RateProcessVI::GetMeanRateVI()  {
          
        double total = 0;
        double* ratealphavi = new double[GetNsite()];
        double* ratebetavi = new double[GetNsite()];
        for (int i=0; i<GetNsite(); i++)  {
                ratealphavi[i] = 0.0;
                ratebetavi[i] = 0.0;
        } 
        GlobalGetRateAlphaVI();
        GlobalGetRateBetaVI();
        for (int i=0; i<GetNsite(); i++)  {
                ratealphavi[i] += RateAlphaVI[i];
                ratebetavi[i] += RateBetaVI[i];
                total += ratealphavi[i]/ratebetavi[i]; 
        }
        return total;
}*/





















































