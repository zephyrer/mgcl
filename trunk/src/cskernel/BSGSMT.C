/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/blcpb.h"
#include "cskernel/SRSmooth2.h"

// THIS IS SURFACE CONSTRUCTION VERSION OF BLGSMT, I.E., 
// RETURN RCOEF AS RCOEF(1,I,.). 
//
// CONSTRUCTS THE CUBIC SMOOTHING SPLINE F TO GIVEN DATA (X(I),Y(I,.)), 
//  I=1,...,NPOINT, WHICH HAS AS SMALL A SECOND DERIVATIVE AS POSSIBLE 
//  WHILE 
//   S(F) 
//    =SUM( ((Y(I,.)-F(X(I)))/DY(I))**2 , I=1,...,NPOINT ) .LE. S, 
//    WHERE S=D*NPOINT 
// ******  I N P U T  ****** 
//  NCD       SPACE DIMENSION OF INPUT DATA ORDINATES Y. 
//  NPOINT.....NUMBER OF DATA POINTS,  A S S U M E D  .GT. 1 
//  X(1),...,X(NPOINT)   DATA POINT  ABSCISSAE,  ASSUMED TO BE STRICTLY 
//        INCREASING . 
//  Y(I,M),...,Y(NPOINT,M)   CORRESPONDING DATA ORDINATES OF NCD SPACE 
//        DIMENSION  (1<=M<=NCD). 
//  DY(I),...,DY(NPOINT)   ESTIMATE OF UNCERTAINTY IN DATA. IF NEGATIVE, 
//        ASSUMED INTERPOLATION REQUIRED. 
//  D.....UPPER BOUND ON THE DISCRETE WEIGHTED MEAN SQUARE DISTANCE OF 
//        THE APPROXIMATION  F  FROM THE DATA . 
//        S=D*NPOINT IS USED AS SUM OF THE DISTANCES 
// ******  W O R K  A R R A Y S  ***** 
//  V.....OF SIZE (NPOINT,7) 
//  A4....OF SIZE (NPOINT,4) 
//  U....OF SIZE (NPOINT,NCD) 
//  PWORK.....OF SIZE (IRC2,NCD) 
// *****  O U T P U T  ***** 
//  N,T(N+4),RCOEF(1,I,NCD)....B-REP OF ORDER 4 OF OBTAINED LINE. 
//        N: B-REP DIMENSION       T(.): KNOT SEQUENCE 
//        RCOEF(1,I,M): B-COEF OF NCD SPACE DIMENSION 
//            1<=I<=N, 1<=M<=NCD 
// ***** NOTE ***** 
//  N=NPOINT+2, AND SO IRC2 MUST BE GREATER THAN OR EQUAL TO (NPOINT+2) 
//  NCD IS ASSUMED LESS OR EQUAL TO 4. 
//
// <<<<< THIS IS STEFFENSEN'S ALGORITHM VERSION  >>>>> 
// STEFFENSEN'S ALGORITHM IS APPLIED TO FIND CORRECT P. 
// ORIGINAL EQUATION TO SOLVE IS: 
//         36*(1-P)**2*DQU - S=0    . 
// THE ALGORITHM IS APPLIED TO THE FOLLOWING VARIATION: 
//         P = 1.- SQRT( (S/36)/DQU  ), 
// WHERE DQU IS SQUARE OF 2-NORM OF D*Q*U . 
void bsgsmt_(int ncd, int npoint,const double *x, 
	const double *y,const double *dy, double d, int iy, int irc1,
	int irc2, double *v, double *a4, double *u, 
	double *pwork, int *n, double *t, double *rcoef)
{
	int i,j,itg,npm1;
	double p,six1mp,sixp,dyi;

	p=SRSmooth2(ncd,npoint,x,y,dy,d,iy,irc2,v,a4,u,pwork);

    // Parameter adjustments 
    u -= npoint+1;
    a4 -= 5;
    v -= npoint+1;
    --dy;
    y -= iy+1;
    rcoef -= irc1*(irc2+1)+1;
    pwork -= irc2+1;;

// COMPUTE POL.COEFFICIENTS FROM  Q*U (IN POUT(.,.) AND A4(1,.)), 
// AND GET B-COEF'S IN RCOEF(.,.) 

    six1mp = (1.f - p) * 6.f;
    sixp = p * 6.f;
    npm1 = npoint - 1;
    itg = 1;
    for(j=1; j<=ncd; ++j){
		for(i=1; i<=npoint; ++i){
			// Computing 2nd power 
		    dyi = dy[i];
		    a4[(i<<2)+1]=y[i+j*iy]-six1mp*(dyi*dyi)*pwork[i+j*irc2];
		}
		for(i=1; i<=npoint; ++i)
			a4[(i<<2)+3]=u[i+j*npoint]*sixp;
		for(i=1; i<=npm1; ++i){
		    a4[(i<<2)+4] = (a4[((i+1)<<2)+3] - a4[(i<<2)+3]) / v[i+(npoint<<2)];
			a4[(i<<2)+2] = (a4[((i+1)<<2)+1] - a4[(i<<2)+1]) / v[i+(npoint<<2)]
				- (a4[(i<<2)+3]+a4[(i<<2)+4]/3.f*v[i+(npoint<<2)])/2.f*v[i+(npoint<<2)];
		}
		blcpb_(x,&a4[5],npm1,4,1,1,npoint,itg,t,&pwork[j*irc2+1],n);
		for(i=1; i<=*n; ++i)
			rcoef[(i+j*irc2)*irc1+1] = pwork[i+j*irc2];
		itg=2;
    }
}
