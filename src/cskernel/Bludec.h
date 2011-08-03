/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUDEC DECREASES B-REP DIMENSION ,APPPROXIMATING THE OLD B-REP. 
// *** INPUT * 
//     ISM   =1 This call of BUDEC is the first call as new data. 
//           =2 CALLED WITH SAME DATA as previous call. 
//              (Only RCOEF1 can be different from the previous call) 
//     K, N1,T1(N1+K),RCOEF1(IRC1,M),IRC1,M 
//            Original B-Spline. 
//     NDEC   : NUMBER OF KNOTS TO DECREASE 
//     IRC2   : ROW DIMENSION OF RCOEF2 
//     WORK1,WORK2 : WHEN ISM=2, PREVIOUS DATA MUST BE INPUT 
// *** OUTPUT * 
//     N2, T2(N2+K), RCOEF2(IRC2,M) 
//             Updated B-Spline. 
//     IFLAG  : = 1 GOOD 
//            : <> 1 ERROR 
// *** WORK   * 
//     WORK1(N2,2K-1) 
//     WORK2(N2) 
// *** NOTE *   . MAXIMUM MULTIPLICITY OF INNER KNOT IS K-1 
void bludec_(int ism, int k, int n1, 
	const double *t1, const double *rcoef1, int irc1, int m, 
	int ndec, int irc2, double *work1, double *work2, 
	int *n2, double *t2, double *rcoef2, int *iflag);
