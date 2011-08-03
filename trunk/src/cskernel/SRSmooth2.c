/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/blgsm1.h"
#include "cskernel/blgsm2.h"
#include "cskernel/blcpb.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//A dedicated subroutine of SRSmooth(Shoengberg and Reinch sommothing function).
//Finds p(the function's return value) and return u[.] and qu[.][.], for
//SRSmooth.
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
//  QTY....OF SIZE (NPOINT,NCD) 
// *****  O U T P U T  ***** 
//  V(NPOINT,7)...PUT DELX = X(.+1) - X(.)  INTO  V(.,4), 
//       PUT THE THREE BANDS OF  Q-TRANSP*D  INTO  V(.,1-3), AND 
//       PUT THE THREE BANDS OF (D*Q)-TRANSP*(D*Q)  AT AND ABOVE THE DIAGONAL
//       INTO  V(.,5-7).
//  U(NPOINT,NCD) 
//  QU..OF SIZE (IQU,NCD), IQU MUST BE GREATER OR EQUAL TO NPOINT
// ***** NOTE ***** 
//  NCD IS ASSUMED LESS OR EQUAL TO 4. 
//
//<<<<< THIS IS STEFFENSEN'S ALGORITHM VERSION OF SMOOTH OF C.DEBOOR >>>>>
// STEFFENSEN'S ALGORITHM IS APPLIED TO FIND CORRECT P. 
// ORIGINAL EQUATION TO SOLVE IS: 
//         36*(1-P)**2*DQU - S=0    . 
// THE ALGORITHM IS APPLIED TO THE FOLLOWING VARIATION: 
//         P = 1.- SQRT( (S/36)/DQU  ), 
// WHERE DQU IS SQUARE OF 2-NORM OF D*Q*U . 
double SRSmooth2(int ncd, int npoint, const double *x, 
	const double *y,const double *dy, double d, int iy, int iqu,
	double *v, double *qty, double *u, double *qu
){
    static const double dpmin = .001;
    static const double divmin = 1e-9;

    // Builtin functions 
    double sqrt(double);

    // Local variables 
    double dbyd, p1mp0;
    int i;
    double p;
    double s, p0, p1, p2,p0old;
    double dp, s36, div;
    double dqu;
    int icount2;
	double ivalues[5]={.5,.25,.75,0.,1.};
	int set0=0, set1=0;

    // Function Body 
    s36=s=0.f;
	if(d>0.f){
		dbyd = d;
		dbyd*=dbyd;
		for(i=0; i<npoint; ++i)
			s+=dbyd/(dy[i] * dy[i]);
		s36 = s/36.f;
	}
    blgsm1_(ncd,x,dy,y,npoint,iy,v,qty);

    if(s>0.f) p=0.f;
    else p = 1.f;
    blgsm2_(p,v,qty,npoint,dy,ncd,iqu,u,qu,&dqu);
    if(s==0.f || dqu<=(s36*1.01f))
		return p;

// ==== STEFFENSEN'S ITERATION STARTS HERE (INITIAL VAL = 0.5) === 
   icount2=0;
   p0 =ivalues[icount2]; p0old=-1.;
L20:
	if(p0==p0old){
	//When  iteration did not converge, we try another initial values, .25 and .75.
		icount2++;
		if(icount2>=5){
			goto L50;
		}
		p0 =ivalues[icount2];set0=set1=0;
	}
L30:
	p0old=p0;
	while(1){
		blgsm2_(p0, v, qty, npoint, dy, ncd, iqu, u,qu, &dqu);
		p1 = 1.f - sqrt(s36/dqu);
		if(p1>0.f)
			break;
	    p0 *= .5f;
	}
 	while(1){
		blgsm2_(p1, v, qty, npoint, dy, ncd, iqu, u,qu,&dqu);
		p2 = 1.f - sqrt(s36 / dqu);
		if(p2>0.f)
			break;
	    p1 *= .5f;
	}

    p1mp0 = p1 - p0;
    div = p2 - p1 - p1mp0;
    if(fabs(div) <= divmin){
		return p2;
	}

    dp = p1mp0 * p1mp0 / div;
    p0 -= dp;
    if(p0<0.f){
		if(set0){
			icount2++;
			if(icount2>=5){	p=p0;goto L50;}
			p0 =ivalues[icount2];set0=0;
			goto L30;
		}
		p0 = 0.f;set0=1;
	}
    if(p0>1.f){
		if(set1){
			icount2++;
			if(icount2>=5){	p=p0;goto L50;}
			p0 =ivalues[icount2];set1=0;
			goto L30;
		}
		p0 = 1.f;set1=1;
	}
    if(fabs(dp)>dpmin)
		goto L20;

// CORRECT VALUE OF P HAS BEEN FOUND. 
L50:
    p = p0;
    blgsm2_(p,v,qty,npoint,dy,ncd,iqu,u,qu,&dqu);
// ===== END OF ITERATION. 
	return p;
}
