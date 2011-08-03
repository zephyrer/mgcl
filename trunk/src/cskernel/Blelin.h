/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLELIN         BLELIN TO GENERATE INPUT DATA OF BLG4SQ FROM B-REP. 
// SUBROUTINE TO EVALUATE B-SPLINE(T(N+K),RCOEF(IRC,NCD)) 
// AT EACH DATA-POINT TAU(I),I=1,--,NT , AND CREATE VAL(NT,NCD) 
// *** INPUT * 
//  K,N,T(N+K),RCOEF(IRC,NCD),IRC,NCD....LINE B-REP TO EVALUATE 
//  IBCBEG.....INDICATES WHAT DATA BE GENERATED AT START 
//          =1 1ST DERIV    =2 2ND DERIV  =3 NO BC 
//  IBCEND.....INDICATES WHAT DATA BE GENERATED AT END 
//          =1 1ST DERIV    =2 2ND DERIV  =3 NO BC 
//  NT,TAU(NT)..DATA-POINTS 
//  IV1,IV2....1ST AND 2ND ARRAY LENGTH OF THE VARIABLE VAL. 
// *** OUTPUT * 
//  VAL(IV1,IV2,NCD).....VALUES EVALUATED 
// *** NOTE * 
//  EVALUATED VALUES ARE STORED IN ROW-WISE IN VAL, I.E. 
//         VAL(1,I,.) I=1....NT. 
void blelin_(int k, int n, const double *t, 
	const double *rcoef, int irc, int ncd, int ibcbeg, int ibcend,
	int nt,const double *tau, int iv1, int iv2, double *val
);
