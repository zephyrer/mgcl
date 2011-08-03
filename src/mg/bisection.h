/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGBISECTION_HH_
#define _MGBISECTION_HH_

///A virtual super class to solve non-linear equations by the bicection method.
///Let f(t) be a function of one double parameter, and a solution of f(t) is known
///to exist between [ts,te]. Then MGBisect::solve() gets the solution, given the initial
///candidate t, the initial span, and the tolerance to halt the iteration.
///To use MGBisection, define a subclass of MGBisection, and prepare member functions
///set_initial_t() and compare_replace(), and use solve() method to get the solution.
class MGBisection{
public:

MGBisection(
	double ts,	///<parameter range from ts
	double te	///<to te.
);

//Virtual Destructor
virtual ~MGBisection(){;};

///Prepare by setting the initial value
virtual void set_initial_t(double t)=0;

///compare with the previous function value(the initial value is set
///by set_initial_t) and replace t with the previous one if necessary.
///The function's return value is the new parameter value.
virtual double compare_replace(
	double t,		///<parameter value to compare at.
	bool& replaced	///<true will be returned if the t is the new solution candidate value.
)=0;

///Solve the equation with bisection method.
double solve(
	double t,	///<The initial parameter value.
	double span_initial,///<The initial parameter span length to increment or decrement,
				///<span*.5 will be the 1st try of the iteration.
	double tolerance,///<The tolerance to halt the bisection iteration.
	int& nrepition	///<iterated number will be returned.
);

private:
	double m_ts, m_te;///<the curve's parameter range from m_ts to m_te;
};

/** @} */ // end of ALGORITHM group

#endif
