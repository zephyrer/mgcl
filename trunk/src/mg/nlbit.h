/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGNLBIT_HH_
#define _MGNLBIT_HH_
/** @addtogroup ALGORITHM
 *  @{
 */

#include <math.h>

/** 
 *	@brief Compute the solution fn(x)=0.
 *  mgNlbit can be applied when known there exists a solution between (xl, xr).
 * SOLUTION OF A NON-LINEAR EQUATION BY BISECTION METHOD AND REGULA FALSI METHOD.
 *
 *	@retval the solution x obtained.
 */
template<class func>
double mgNlbit(
	func& fn,	///<Function object to evaluate the fucntion,
				///<must return double value.
	double xl,	///<LEFT HAND SIDE VALUE OF THE SOLUTION.
	double xr,	///<RIGHT HAND SIDE VALUE OF THE SOLUTION  . 
	double eps,	///<tolerance allowed for the convergence in world coordnate.
	int itr,	///<MAXIMUM NUMBER OF REPITITION.
	int& ier	///<Error code. =0:solution successfully obtained,
				///<     =1:mgNlbit did not converg, and solution not obtained.
){
    ier = 0;
	double epsHalf = eps*.5;
    double x = xl;
    double fl = fn(x);
    if (fabs(fl)<=eps) return x;
    x = xr;
    double fr = fn(x);
    if(fabs(fr)<=eps) return x;
    if(fr * fl > 0.) {
		x = (xl + xr) * .5;
		ier = 2;
		return x;
    }

    for(int n=1; n<=itr; ++n){
		x=(xl+xr)*.5;
		double f = fn(x);
		if(fabs(f)<=eps)
			return x;

		if(f*fl < 0.){
		    xr = x;
			fr = f;
		}else{
		    xl = x;
			fl = f;
		}
		double df=fl-fr;
		if(fabs(df)<=epsHalf){
			x=(xr+xl)*.5;
			return x;
		}
		double dx = (xr-xl)*fl/df;
		x = xl+dx;
		f = fn(x);
		if (fabs(f)<=eps)
			return x;

		if(fl*f<0.){
		    xr = x;
			fr = f;
		}else{
		    xl = x;
			fl = f;
		}
    }
    ier = 1;
    return x;
}

/** @} */ // end of ALGORITHM group
#endif
