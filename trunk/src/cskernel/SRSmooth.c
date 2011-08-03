/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/SRSmooth2.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Shoengberg and Reinch sommothing function.
//Instead of computing B-Spline, SRSmooth will compute only smoothed points
//in pout.
//
//  CONSTRUCTS THE CUBIC SMOOTHING SPLINE F TO GIVEN DATA (X(I),Y(I,.)), 
//  I=1,...,NPOINT, WHICH HAS AS SMALL A SECOND DERIVATIVE AS POSSIBLE 
//  WHILE 
//  S(F) = SUM( ((Y(I,.)-F(X(I)))/DY(I))**2 , I=1,...,NPOINT ) .LE. S, 
//  WHERE S=D*NPOINT 
// ******  I N P U T  ****** 
//  NCD       SPACE DIMENSION OF INPUT DATA ORDINATES Y. 
//  NPOINT.....NUMBER OF DATA POINTS,  A S S U M E D  .GT. 1 
//  X(1),...,X(NPOINT)   DATA POINT  ABSCISSAE,  ASSUMED TO BE STRICTLY 
//        INCREASING . 
//  Y(1,J),...,Y(NPOINT,J)   CORRESPONDING DATA ORDINATES OF NCD SPACE 
//        DIMENSION  (1<=J<=NCD). 
//        Y(IY,NCD) is the array length.
//  DY(1),...,DY(NPOINT)   ESTIMATE OF UNCERTAINTY IN DATA. IF NEGATIVE, 
//        ASSUMED INTERPOLATION REQUIRED. 
//  D.....UPPER BOUND ON THE DISCRETE WEIGHTED MEAN SQUARE DISTANCE OF 
//        THE APPROXIMATION  F  FROM THE DATA . 
//        S=D*NPOINT IS USED AS SUM OF THE DISTANCES 
// ******  W O R K  A R R A Y S  ***** 
//  V.....OF SIZE (NPOINT,7) 
//  A4....OF SIZE (NPOINT,NCD) 
//  U....OF SIZE (NPOINT,NCD) 
//  WORK..OF SIZE (NPOINT,NCD)
// *****  O U T P U T  ***** 
//  POUT(IPO,NCD)....CONTAINS THE SEQUENCE OF SMOOTHED ORDINATES . 
//        THE LENGTH OF POUT IS NPOINT.
// ***** NOTE ***** 
//  NCD IS ASSUMED LESS OR EQUAL TO 4. 
//  Y(.,.) AND POUT(.,.) MAY BE THE SAME AREA. 
//
//<<<<< THIS IS STEFFENSEN'S ALGORITHM VERSION OF SMOOTH OF C.DEBOOR >>>>>
// STEFFENSEN'S ALGORITHM IS APPLIED TO FIND CORRECT P. 
// ORIGINAL EQUATION TO SOLVE IS: 
//         36*(1-P)**2*DQU - S=0    . 
// THE ALGORITHM IS APPLIED TO THE FOLLOWING VARIATION: 
//         P = 1.- SQRT( (S/36)/DQU  ), 
// WHERE DQU IS SQUARE OF 2-NORM OF D*Q*U . 
void SRSmooth(int ncd, int npoint, const double *x, 
	const double *y,const double *dy, double d, int iy, int ipo,
	double *v, double *a4, double *u, double *work,	double *pout
){
	double six1mp,dyi;
	int i,j;
	double p=SRSmooth2(ncd,npoint,x,y,dy,d,iy,npoint,v,a4,u,work);

    pout -= ipo + 1;
    --dy;
    y -= iy + 1;
    work -= npoint + 1;
// COMPUTE SMOOTHED POINTS IN POUT[] 
    six1mp = (1.f - p) * 6.f;
    for(j=1; j<=ncd; ++j){
		for(i=1; i<=npoint; ++i){
		    dyi = dy[i];
			pout[i+j*ipo]=y[i+j*iy]-six1mp*(dyi*dyi)*work[i+j*npoint];
		}
    }
}
