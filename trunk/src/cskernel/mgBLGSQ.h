/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  mgblgsq IS EASY-TO-USE VERSION OF blg4sp2_, I.E. mgblgsq GENERATES 
//  KNOT VECTOR USING blg4s1_, THEN B-COEF'S USING blg4sp2_. 
// *** INPUT * 
//         k......order of the b-spline.
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
//         T(N+k) : KNOT VECTOR 
//         RCOEF(IRC,NCD) : B-SPLINE 
//         IFLAG  : =1 NORMAL END 
//                  =2 ABNORMAL 
// *** WORK * 
//         WORK(N,k*2+1) : WORK AREA FOR SUBROUTINE blg4sp2_ 

// ****    CREATE KNOT-VECTOR T(N+4) FROM TAU(N) 
int mgblgsq(int k, int ibcbeg, int ibcend, double *tau,
	double *val, int iv, int n, int ncd, int irc,
	double *work, double *t, double *rcoef, int* iflag);
