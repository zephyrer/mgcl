/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/bisection.h"

//A virtual super class to solve non-linear equations by bicection methos.
MGBisection::MGBisection(
	double ts,	//parameter range from ts to te.
	double te
):m_ts(ts), m_te(te){
}

//Compute the fn(t)'s parameter value that is the maxima 
//Function's return value will be the solution obtained.
double MGBisection::solve(
	double t,	//The initial parameter value.
	double span_initial,//The initial parameter span length to increment or decrement.
				//span*.5 will be the 1st try of the iteration.
	double tolerance,//The tolerance to halt the bisection iteration.
	int& nrepition	//iterated number will be returned.
){
	assert(tolerance>0.);

	double span=span_initial;
	nrepition=0;
	set_initial_t(t);
	bool replaced;
	while(span>=tolerance){
		span*=.5; nrepition++;
	    double tr = t+span;
		if(tr<m_te){
			t=compare_replace(tr, replaced);
			if(replaced)
				continue;
		}
	    double tl = t-span;
		if(tl>m_ts)
			t=compare_replace(tl, replaced);
	}
    return t;
}
