/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  REAL FUNCTION TO GET INTGERAL OF F(T)=X*DY-Y*DX. THIS INTEGRAL IS 
// FOR AREA COMPUTATION, I.E. AREA BOUNDED BY B-REP. 
// *** INPUT * 
// K,N,T(N+K),RCOEF(IRC,2).....B-REP OF ORDER K AND 2 SPACE DIMENSION. 
// T1,T2.....PARAMETER RANGE OF THE INTEGRATION, FROM T1 TO T2. 
// *** WORK * 
// TW(MM+2*K-2),WK1(MM),WK2(MM,4*K-4)....WORK AREA OF EACH LENGTH. HERE 
//         MM=N IF K=2, OTHERWISE MM=(N-K)*K+4*K-4. 
// *** OUTPUT * 
// BLA2D....THE INTEGRAL OF THE FUNCTION (X*DY-Y*DX), HERE X=RCOEF(.,1), 
//         Y=RCOEF(.,2), DX=DX/DT, DY=DY/DT. 
double bla2d_(
	int k, int n,const double *t,const double *rcoef, 
	int irc, double t1, double t2, double *tw, 
	double *wk1, double *wk2);
