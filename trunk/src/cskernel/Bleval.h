/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLEVAL EVALUATES THE JDERIV-TH DERIVATIVE(S), GIVEN B-COEF'S 
// RCOEF IN COLUMN-WISE OR ROW-WISE. B-COEFFICIENTS RCOEF MAY 
// BE MULTIPLE OF THE SAME KNOT CONFIGURATION. 
// ***INPUT* 
//    K,N,T(N+K),RCOEF(IRC,.),IRC,NCD,KAR....DESCRIBE B-REP. 
//           ORDER, B-REP DIMENSION, KNOT VECTOR, B-COEF, AND SPACE 
//           DIMENSION. B-COEF MAY BE STORED IN COLUMN-WISE OR IN ROW- 
//           WISE ACCORDING TO KAR. I.E. 
//           KAR= 1 :  RCOEF(I,J)   1<=I<=N  (J=1...NCD) (COLUMN-WISE)
//              <>1 :RCOEF(J,I)   1<=I<=N  (J=1...NCD)  (ROW-WISE) 
//    X,JDERIV....PARAMETER VALUE AT WHICH DERIVATIVE TO EVALUATE AND 
//           ORDER OF DERIVATIVE. THE DERIVATIVE ORDER JDERIV MUST BE 
//           NON-NEGATIVE, MAY BE ZERO. 
//    JCONT.......INDICATES WHICH KIND OF CONTINUITY BE REQUIRED. I.E. 
//             WHEN JCONT= 1 : RIGHT-CONTINUOUS AT PARAMETER X 
//                       <>1 : LEFT-CONTINUOUS. 
//    ISM....INDICATES WHETHER OR NOT BLEVAL SHOULD START THE PROCESS 
//          FROM THE EVALUATION OF B-SPLINE (ISM=1) OR THE EVALUATION 
//          CAN BE OMITTED BECAUSE THIS CALL IS THE SECOND CALL OF THE 
//          SAME KNOT AND THE SAME PARAMETER VALUE X,AND JDERIV (ISM<>1). 
// ***OUTPUT* 
//    P(NCD)         EVALUATED DERIVATIVE(S) 
// ***NOTE* 
//  . BLEVAL IS MAINLY FOR MULTIPLE B-COEF'S OF THE SAME KNOT VECTOR, 
//    B-COEFS OF ROW-WISE, OR LEFT-CONTINUOUS EVALUATION. IF THE OBJECT 
//    B-REP IS SIMPLE ENOUGH, BLE IS EASY TO USE AND HAVE BETTER 
//    PERFORMANCE. SEE BLE. 
int bleval_(
	int k, int n,const double *t,const double *rcoef, int irc, int ncd, int kar, 
	double x, int jderiv, int jcont, int ism, 
	double *p
);
