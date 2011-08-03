/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//Get the relative zero tolerance.
double bzrzro_();

//Get the absolute zero tolerance;
double bzamin_();

//Get the machine zero tolerance;
double bzmzro_();

//Get the maximum ratio of adjacent knot vector spans.
double bkmax_();

//Print B-Spline data.
void bzprintBspl(int k, int n, const double *t,const double *rcoef, int irc, int ncd);
