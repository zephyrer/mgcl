/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  BLG4SQ IS EASY-TO-USE VERSION OF BLG4SP, I.E. BLG4SQ GENERATES 
//  KNOT VECTOR USING BLG4S1, THEN B-COEF'S USING BLG4SP. 
// *** INPUT * 
//         IBCBEG, IBCEND......BOUNDARY COND OF BEGINNING AND ENDING 
//                POINT, EACH. 
//                =1 1ST DERIV PROVIDED,    =2 2ND DERIV 
//                =3 NO BOUNDARY COND. 
//                =4 BOTH 1ST AND 2ND DDERIVATIVES PROVIDED. 
//         TAU(N) : DATA-POINTS 
//         VAL(IV,NCD) : KOSIN-OFFSET-DATA OF M-ORDER 
//         IV     : COLUMN-LENGTH OF VAL 
//         N      : DATA-NO. OF VAL      ( IV >= N ) 
//         NCD    : ORDER OF COORDINATES 
//         IRC    : COLUMN-LENGTH OF RCOEF 
// *** OUTPUT * 
//         T(N+4) : KNOT VECTOR 
//         RCOEF(IRC,NCD) : B-SPLINE 
//         Function's return value: =1 NORMAL END, =2 ABNORMAL
// *** WORK * 
//         WORK(N,9) : WORK AREA FOR SUBROUTINE BLG4SP 
int blg4sq_(int ibcbeg, int ibcend, const double *tau, const double *val,
	int iv, int n, int ncd, int irc,double *work, double *t, double *rcoef);
