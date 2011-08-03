/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUMOV MOVES A LINE B-REP AS SPECIFIED IN KT. 
// *** INPUT * 
//       K,N,T(N+K),RCOEFI(IRCI,NCD),IRCI,NCD....DESCRIBE B-REP INPUT 
//             K:ORDER,     NCD:SPACE DIMENSIONN,   N:B-REP DIMENSION 
//                T(N+K):KNOT VECTOR,  RCOEFI(.,.): B-COEFFICIENTS 
//       TC...IS THE CENTER OF TRANSLATION, PARAMETER VALUE OF THE B-REP
//       P(NCD).....IS THE NEW POINT COORDINATE, 
//                THE POINT F(TC) IS MOVED TO P(.), WHERE F IS B-REP 
//       KT.....INDICATES WHAT KIND OF MOVE IS BEING PERFORMED. 
//          =1 : TWO END POINTS FIXED AND CENTER OF MOVE BETWEEN THEM. 
//          =2 : ONE POINT,TFI(1),FIXED. THE OTHER POINT IS FREE 
//          =3 : TWO POINTS,TFI(.),FIXED AND CENTER OF MOVE BETWN THEM. 
//          =4 : NO FIXED POINT SPECIFIED, MINIMUM MOVE 
//          =5 : PARALLEL TRANSLATION OF THE LINE 
//       TFI(2)....SPECIFIES PARAMETER VALUES OF B-REP TO FIX AT. 
//          WHEN KT=2, TFI(2) IS DUMMY. 
//          WHEN KT=1, 4 AND 5, TFI(.) ARE DUMMY. 
//       IRC.....ROW DIMENSIONS OF THE VARIABLES RCOEF 
// *** OUTPUT * 
//       RCOEF(IRC,NCD).....ARE NEW TRANSLATED B-COEF'S. 
//          RCOEF MAY BE THE SAME AREA AS RCOEFI 
//       TFO(2).....ARE PARAMETER VALUES OF THE B-REP, GIVES TWO 
//          BOUNDARY POINTS OF THE MOVE. I.E., MOVE IS PERFORMED IN 
//          TFO(1) < T < TFO(2). 
// ***WORK* 
//       WRATIO(N)...WORK AREA OF LENGTH N. 
void blumov_(int k, int n, const double *t, 
	const double *rcoefi, int irci, int ncd, double tc, 
	const double *p, int kt,const double *tfi, int irc, double *wratio,
	double *rcoef, double *tfo);
