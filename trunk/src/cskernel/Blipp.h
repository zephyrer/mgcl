/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLIPP COMPUTES INTERSECTION PARAMETER VALUE X(I) OF 1-D B-REP, S.T. 
//  F = G(X(I)) FOR G(.):B-REP FUNCTION, AND X(I) > TS FOR 1<=I<=MX. 
// *** INPUT * 
//     K,N,T(N+K),RCOEF(N)......BREP OF ONE SPACE DIMENSION. 
//     F......  B-REP FUNCTION VALUE TO GET THE PARAM VALUES. 
//     ERROR.....ERROR ESTIMATE ALLOWED TO GET THE INTERSECTION POINTS. 
//     TS.......INDICATES PARAMETER VALUE AT WHICH TO START THE 
//              COMPUTATION. 
//     MX......PROVIDES LENGTH OF THE VARIABLE X(.), BLIPP STOPS 
//             COMPUTATION WHEN NUM OF INTERSECTION POINTS REACHES MX. 
// *** OUTPUT * 
//     NX,X(NX)..... PARAMETER VALUES OBTAINED BY BLIPP. NX<=MX. 
//                   NX=0 WHEN NO SOLUTION. 
//     IEND.....GIVES THE INF WETHER BLIPP COMPUTES TO THE END OF 
//              PARAMETER. 
//               =0: NOT TO THE END BECAUSE MX IS TOO SMALL. 
//               =1: COMPUTE TO THE END. 
// *** WORK * 
//     WORK(.)       WORK ARRAY OF LENGTH 4*K*K+3*K 
void blipp_(int k, int n, const double *t, const double *rcoef,
	double f, double error, double ts, int mx,
	double *work, int *nx, double *x, int *iend);
