/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
);
