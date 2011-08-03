/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//                                                     Y. MIZUNO 
// BLUMI2 IS AN INTERNAL SUBROUTINE OF BLUMIX, 
// BLUMI2 WILL GET PARAM VALUE OF B-REP (N,T,RCOEF) WHOSE FUNCTION VALUE 
// IS F(.). 
// *** INPUT  * 
//     KCOD1,K,N,T(.),RCOEF(IRC,2),IRC....DESCRIBE THE B-REP OF ORDER K. 
//     KCOD2,TAU,F(2)......ARE DATA POINT AND ASSOCIATED FUNCTION VALUE 
//            , KCOD2 SPECIFIES COORDINATE KIND. 
//     TPREV.....IS THE PARAM VALUE USED WHEN MORE THAN TWO INTERSECTION 
//            POINTS, NEAREST AND GREATER THAN TPREV PARAM VALUE IS 
//            EMPLOYED. 
// *** OUTPUT *** 
//     TINT......PARAM VALUE OF THE B-REP OBTAINED. 
// ***WORK * WORK(4*K*K+3K)....WORK OF LENGTH 4*K*K+3K. 
void blumi2_(int kcod1, int k, int n, 
	const double *t, const double *rcoef, int irc, int kcod2, 
	double tau, double *f, double tprev, double error, 
	double *work, double *tint);

// BLUMI4         BLUMI4 FOR BLUMIX 
// BLUMI4 GETS UNIT TANGENT VECTOR FROM TWO TWO-DIMENSIONAL TANGENT 
// VECTORS F1(2) AND F2(2) AND STORE THEM IN RCOEF(1,J) FOR 1<=J<=3. 
//   *** KCOD1 MUST NOT BE 4 *** 
void blumi4_(int kcod1, double *f1, int kcod2, 
	double *f2, double *rcoef, int irc);
