/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUPRT computes a partial b-spline, given a spline and 
// a partial parameter range of the original knot configuration. 
// The new spline is exactly the same as the old one, 
// although the representation is partial. 
// *** INPUT * 
//     K,N1,T1(N1+K),RCOE1(IRC1,NCD),IRC1,NCD... 
//                     Describe the original B-REP. 
//     TSI,TEI.....New parameter range, should be 
//               TSI < TEI.  If TSI <= T1(K) or T1(N1+1) <= TEI, 
//               obtained b-rep is not partial and the same as 
//               the original one, regarding to the start or end 
//               side of the origial line each. 
//     IRC2....ROW DIMENSION OF THE VARIABLE RCOE2 
//     multiple...indicates if knot multiplicity of K is necessary at start
//                and end parameter of T2.
//               =0: unnecessary, !=0: necessary.
// *** OUTPUT * 
//     N2,T2(N2+K),RCOE2(IRC2,NCD)..THE NEW B-COEF OBTAINED. 
// *** WORK * 
//     WORK(K,K)  LENGTH OF K*K 
// *** NOTE * 
//     If TSI >= TEI, N2=0 is set and nothing will be done. 
void bluprt_(int k, int n1, const double *t1, 
	const double *rcoe1, int irc1, int ncd, double tsi, 
	double tei, int irc2, double *work, int *n2, 
	double *t2, double *rcoe2, int multiple);
