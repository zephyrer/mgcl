/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLE EVALUATES THE JDERIV-TH DERIVATIVE(S), GIVEN B-COEF'S 
// RCOEF IN COLUMN-WISE OR ROW-WISE. B-COEFFICIENTS RCOEF MAY 
// BE MULTIPLE OF THE SAME KNOT CONFIGURATION. 
// ***INPUT* 
//    K,N,T(N+K),RCOEF(IRC,.),IRC,NCD....DESCRIBE B-REP. 
//           ORDER, B-REP DIMENSION, KNOT VECTOR, B-COEF, AND SPACE 
//           DIMENSION. 
//             RCOEF(I,J)   1<=I<=N  (J=1...NCD) 
//    X,JDERIV....PARAMETER VALUE AT WHICH DERIVATIVE TO EVALUATE AND 
//           ORDER OF DERIVATIVE. THE DERIVATIVE ORDER JDERIV MUST BE 
//           NON-NEGATIVE, MAY BE ZERO. 
// ***OUTPUT* 
//    P(NCD)         EVALUATED DERIVATIVE(S) 
// ***NOTE* 
//    BLE IS EASY-TO-USE VERSION OF BLEVAL, AND HAVE BETTER PERFORMANCE 
//    THAN BLEVAL. 
void ble_(
	int k,int n,const double *t,const double *rcoef,int irc,int ncd,double x,int jderiv,double *p
);
