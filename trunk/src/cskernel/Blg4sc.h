/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  BLG4SC GENERATES CIRCULAR B-SPLINE, THAT IS, STARTING AND ENDING 
//  POINTS COINCIDE AND THEIR TANGENT'S ARE EQUAL. 
// *** INPUT * 
//         NV,VAL(IV,NCD),IV,NCD.....ORIGINAL DATA OF SPACE DIMENSION 
//                NCD. NV IS NUM OF DATA AS VAL(I,.) 1<=I<=NV. 
//         IRC....ROW-DIMENSION OF RCOEF AS RCOEF(IRC,NCD). 
// *** OUTPUT * 
//         N,T(N+4),RCOEF(IRC,NCD).....B-SPLINE OF ORDER 4 OBTAINED. 
//                THE SPACE DIMENSION IS NCD AND THE B-REP DIMENSION 
//                IS N=NV+2. 
//         IFLAG  : =1 NORMAL END 
//                  <>1 ABNORMAL 
// *** WORK * 
//         TAU(N),WORK(N,9) : WORK AREA FOR SUBROUTINE BLG4SQ 
void blg4sc_(int nv,double *val, int iv, 
	int ncd, int irc, double *tau, double *work, 
	int *n, double *t, double *rcoef, int *iflag
);
