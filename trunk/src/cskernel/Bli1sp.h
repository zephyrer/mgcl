/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// REAL FUNCTION TO GET INTERSECTION POINT OF 1-DIMENSIONAL B-REP. 
// B-REP. IS FIRST CONVERTED INTO PP-REP, THEN NEWTON-RAPHSON METHOD IS
//  USED FOR THE SOLUTION. EVALUATION OF B-REP IS DONE USING PP-REP. 
// *** INPUT *
//   K,N,T(N+K),RCOEF(N).......B-REP FOR INTERSECTION COMPUTATION. 
//   KI........KNOT INDEX OF T S.T. 
//                        BLER(,T(KI),) <= F <= BLER(,T(KI+1),)   , 
//                    OR  BLER(,T(KI),) >= F >= BLER(,T(KI+1),) 
//   F.......THE B-VALUE TO FIND THE ASSOCIATED PARAMETER VALUE BLI1SP, 
//                                  F = BLER(,BLI1SP,). 
//         ERROR               ERROR ESTIMATE OF F, I.E. THE VALUE S.T. 
//  ERROR......ERROR ALLOWED TO COMPUTE INTERSECTION. 
// *** OUTPUT * 
//    BLI1SP........THE PARAMETER VALUE S.T.   F=BLER(,BLI1SP,) 
// *** WORK * 
//   WORK(K*K)         WORK AREA OF K*K 
double bli1sp_(int k, int n, const double *t, const double *rcoef, 
	int ki, double f, double error, double *work);
