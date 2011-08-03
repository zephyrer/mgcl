/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//   B1A COMPUTES THE INTEGRAL OF B-REP (K,N,T,RCOEF) FROM PARAMETER 
//   T(1) TO X. 
//***INPUT*
//  K,N,T(N+K),RCOEF(N)......PROVIDE B-REP(OF ONE SPACE DIMENSION), 
//             FUNCTION TO INTEGRATE. 
//  X.....SPECIFIY PARAMETER RANGE OF THE INTEGRATION, FROM T(1) TO X.
//      X MUST LIE IN THE PARAMETER RANGE OF THE B-REP, I.E. 
//                    T(1) <= X <= T(N+K). 
//***OUTPUT*
//  B1A......INTEGRAL FROM PARAMETER T(1) TO X. 
double b1a_(int k, int n,const double *t,const double *rcoef, double x);
