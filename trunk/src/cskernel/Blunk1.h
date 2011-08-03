/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//                                                 Y.MIZUNO 
//         BLUNK1 WILL GET THE J-TH B-COEFFICIENTS, GIVEN NEW KNOT 
//         CONFIGURATION T(I), USED WHEN ADDING NEW PSEUDO KNOT. 
// ** INPUT * 
//     TAU(N+K)     : THE OLD KNOT VECTOR          N=B-REP DIMENSION OF 
//     RCOEF(N)     : THE OLD B-COEFFICIENTS.        THE OLD. 
//     K            : THE ORDER OF B-REP. 
//     MU           : THE SUBSCRIPT OF TAU S.T. 
//                       TAU(MU) <= T(J) < TAU(MU+1) . 
//                    (  SHOULD NOT  TAU(MU)=TAU(MU+1) ) 
//     T(M+K)       : THE NEW KNOT VECTOR 
//     J            : INDICATES WHICH B-COEF. SHOULD BE EVALUATED AMONG 
//                    M B-COEF. OF THE NEW B-REP. 
// ** OUTPUT * 
//     BLUNK1 = J-TH B-COEFFICIENTS OF THE NEW B-REP. 
// ** WORK * 
//     WORK(K*K)     WORK ARRAY OF LENGTH K*K 
// ** NOTE * 
//     BLUNK1 EMPLOYS THE OSLO ALGORITHM OF COHEN,ET AL. 
double blunk1_(
	const double *tau, const double *rcoef, int k, int mu, const double *t, int j,
	double *work
);
