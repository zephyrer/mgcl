/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//         BLIPP3 WILL OBTAIN THE SOLUTION OF F=G(X) -- G:B-REP-- BY 
//         CONVERTING EACH SPAN INTO BEZIER CURVE AND SUBDIVIDING THE 
//         BEZIER CURVE. 
//         THE SOLUTION X (PARAM VALUE) ARE MORE THAN ONE,IN GENERAL. 
//         OSLO ALGORITHM IS EMPLOYED FOR SUBDIVISION. 
// *** INPUT * 
//     K,N,T(N+K),RCOEF(N).....PROVIDE B-REP. 
//     KI            INDEX OF RCOEF WHERE BLIPP3 TO START GETTING THE 
//                   SOLUTION. THE SOLUTION X IS >= T(KI+1) 
//     F             FUNCTION VALUE TO GET THE SOLUTION OF F=G(X) 
//     ERROR         ERROR ESTIMATE S.T. ABS(F-G(X)) <= ERROR 
//     ERRORT        ERROR ESTIMATE OF KNOT VECTOR. 
//     MX............MAX LENGTH OF TH EVARIABLE X(.). 
//     TS...........LOWER PARAMETER LIMIT, I.E. OBTAINED X(.)>TS 
//                  ( X(.) SUCH THAT X(.)<=TSIS DISCARDED. ) 
//     NX            NUMBER OF THE SOLUTIONS OBTAINED BEFORE THE CALL 
//                   OF BLIPP3. NUMBER OF NEW SOLUTIONS WILL BE ADDED 
//                   TO NX. 
// *** OUTPUT * 
//     NX,X(NX)..... THE SOLUTIONS ARE SET INTO X(I) AND NX UPDATED. 
//     KI............INDEX OF RCOEF IS RETURNED S.T. 
//                    RCOEF(KI-I) > F FOR I=0,--,K-2    , 
//                OR  RCOEF(KI-I) < F FOR I=0,--,K-2. 
//                   WHEN KI=N, COMPUTATION REACHED TO THE END. 
//     IEND..........INDICATES IF BLIPP3 PERFORMED FULL COMPUTATION 
//                  OR NOT. IEND=0 NOT, BECAUSE MX IS TOO SMALL. 
//                          IEND=1 PERFORMED FULL COMPUTATION 
// *** WORK * 
//     PSTAK(K,K+1)  WORK ARRAY OF LENGTH K*(K+1). 
//                   THE FIRST K*2 WORDS FOR TEMPORAL BEZIER CURVE. 
//                   THE LAST K*(K-1) WORDS FOR BEZIER CURVE STACK AREA. 
//     TSTAK(K,2,K+1) USED FOR KNOT VECTOR AREA OF THE PSTAK. 
//                   THE LAST K*(K-1) WORDS FOR BEZIER CURVE STACK AREA. 
//     WORK3(K*K) WORK ARRAY OF LENGTH K*K, FOR BLUNK1 AND BLI1SP. 
//                     USED FOR RCBNK1,BLI1SP. 
void blipp3_(int kin, int n,const double *t, 
	const double *rcoef, int *ki, double f, double error, 
	double errort, int mx, double ts, double *pstak, 
	double *tstak, double *work3, int *nx, double *x, 
	int *iend);

