/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCurveParameter_HH_
#define _MGCurveParameter_HH_

#include <math.h>
#include "mg/nlbit.h"
#include "mg/Interval.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/** @addtogroup ALGORITHM
 *  @{
 */

///Utility class to compute a curve parameter defined by f(t)=0.
///f(t) is defined by "virtual double operator()=0;".
///To use MGCurveParameter class, define subclass of MGCurveParameter,
///and define operator().
///Then use getCurveParameter member function to compute the solution of f(t)=0.
///Before using getCurveParameter(), the tolerance and the incremental value must be
///set. These default values provided in the constructor are usually not valid ones.
class MGCurveParameter{

public:

///constructor. error and delta must be set if not specified in this
///constructor. Following default values are usually not effective.
MGCurveParameter(
	const MGInterval& prange,///<parameter range for the computation.
	double error=.001,	///<tolerance to compute the parameter value(for f(t)).
	double delta=5.	///<incremental value for the curve parameter.
):m_prange(prange), m_error(error), m_delta(delta){;};

///return the function value f(t) to solve f(t)=0.
virtual double operator()(double t)const=0;

///Get the paramer value of f(t)=0 that is defined by operator()(double)=0.
///Function's return code is:
///0: when the solution is successfully obtaine in t,
///1: There is no solution.
///-2:system error(usually this does not occur. If occured, some bugs are included.)
int getCurveParameter(
	double& t	///<input the guess parameter, the exact solution
		///<will be returned when function's return value is 0.
);

void set_delta(double delta){m_delta=delta;};
void set_error(double error){m_error=error;};

private:
	const MGInterval m_prange;	///<parameter range for the computation.
	double m_error;				///<tolerance to compute the parameter value.
	double m_delta;				///<incremental value for the curve parameter.

///Get the paramer value of f(t)=0 that is defined by operator()(double)=0.
///Function's return code is:
///0: when the solution is successfully obtaine in t,
///1: solution not obtained for the guess parameter t(try with other guess parameter).
///-2:system error(usually this does not occur. If occured, some bugs are included.)
int MGCurveParameter::getCurveParameter2(
	double& t	///<input the guess parameter, the exact solution
		///<will be returned when function's return value is 0.
)const;
};

/** @} */ // end of ALGORITHM group

#endif