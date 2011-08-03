/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//         SUBROUTINE TO PRODUCE THE B-SPLINE COEF.S RCOEF OF THE SPLINE 
//         OF ORDER k WITH KNOTS T(I),I=1,2,--,N+K , 
//         GIVEN DATA POINTS TAU(I) AND ASSOCIATED DATA VAL(I,.). 
//         THE   VAL(I,J) TAKES POSITIONA DATA AT TAU(I) 
//                   FOR TAU(I-1) < TAU(I) < TAU(I+1) , 
//             VAL(I,J) = D(F(TAU(I)), VAL(I+1,J) = F(TAU(I+1)), 
//                   FOR TAU(I-1)<TAU(I)=TAU(I+1)<TAU(I+2), 
//             VAL(I,J) = D(F(TAU(I)), VAL(I+1,J)=F(TAU(I+1)), 
//             VAL(I+2,J)=D(F(TAU(I+2)), 
//                   FOR TAU(I-1)<TAU(I)=TAU(I+1)=TAU(I+2)<TAU(I+3). 
// *** INPUT * 
//         K         ORDER OF THE SPLINE TO MAKE.
//         IFLAG     SPECIFIES BLG4SP2 TO PRODUCE Q,WORK2 FROM B.C. 
//                   AND T(KNOT) -- IFLAG=2 --, OR Q AND WORK2 CAN 
//                   BE UTILIZED BECAUSE THE PREVIOUS CALL OF BLG4SP2
//                   PRODUCED THEM AND THIS CALL IS THE SECOND CALL 
//                   WITH THE SAME KNOT AND BOUNDARY COND 
//                   -- IFLAG=1 --  . 
//         IBCBEG    SPECIFY BOUNDARY COND. AT T=T(K) (IBCBEG), AND 
//         IBCEND    AT T(N+1) (IBCEND) . 
//                =1 : FIRST DERV. PROVIDED IN VAL(2,J) OR VAL(N-1,J) 
//                =2 : 2ND DERIV. PROVIDED IN VAL(2,J) OR VAL(N-1,J) 
//                =3 : NO BOUNDARY COND. 
//                =4 : 1ST AND 2ND DERV. PROVIDED IN (VAL(2,J),VAL(3,J)) 
//                      , OR IN (VAL(N-1,J),VAL(N-2,J)) EACH 
//         TAU(N)    DATA POINTS SEQUENCE OF LENGTH N. 
//                   TAU(I)<=TAU(I+1)  (AT MOST THREE MULTIPLICITIES). 
//         VAL(IV,M) PROVIDES INTERPOLATION DATA AS ABOVE. 
//                   HAS DIMENSION  VAL(IV,M),AND VAL(I,J),I=1,--,N 
//                   PROVIDES DATA. 
//         IV        ROW DIMENSION OF VAL SPECIFIED IN THE CALLING 
//                   PROGRAM. 
//         N         IS THE B-REP. DIMENSION. 
//         M         SPECIFIES NUMBER OF DATA SETS, I.E. NUMBER OF 
//                   CURVES TO BE INTERPOLATED WITH THE SAME KNOTS. 
//         T(N+K)    KNOT VECTOR OF THE B-REP. 
//         IRC       SPECIFIES THE ROW DIMENSION OF RCOEF. 
//         Q(2*K-1,N)    THE SAME CONTENTS OF THE PREVIOUS CALL. 
//                   (NECESSARY ONLY WHEN IFLAG=1) 
//         WORK2(N)  THE SAME CONTENTS OF THE PREVIOUS CALL. 
//                   (NECESSARY ONLY WHEN IFLAG=1) 
// *** WORK * 
//    WORK1(N)   ARRAY OF LENGTH N . 
//    WORK2(N)   ARRAY OF LENGTH N . 
//                   (WHEN IFLAG=1, INPUT) 
// *** OUTPUT * 
//         Q(2*K-1,N)  ARRAY OF SIZE (2K-1)*N CONTAINING THE TRIANGULAR 
//                   FACTORIZATION OF THE COEFFICIENT MATRIX OF THE 
//                   LINEAR SYSTEM.(WHEN IFLAG=2) 
//         RCOEF(IRC,N) B-COEFFICIENTS OF THE INTERPOLANT,PROVIDED IN 
//                   ROW-WISE. 
//         IFLAG     INDICATES SUCCESS (=1) OR FAILURE(=2). 
// *** NOTE * 
//        WHEN IFLAG=1, IBCBEG,IBCEND,TAU, AND,T WILL NOT BE USED, 
//        INSTEAD Q AND WORK2 ARE NECESSARY AS INPUT. 
//        DATA POINTS AND KNOT VECTOR HAVE THE FOLLOWING RESTRICTION 
//      . TAU(1)=T(K) , TAU(N)=T(N+1) 
//      . WHEN TAU(I)=TAU(I+1) (=TAU(I+2)), 
//             T(I+2)=T(I+3) (=T(I+4)) =TAU(I). 
//      . T(1)=...=T(K) ,  T(N+1)=...=T(N+K). 
void blg4sp2_(int k, int *iflag, int ibcbeg, int ibcend,
	 const double *tau, const double *val, int iv, int n, int m,const double *t, int irc,
	 double *work1, double *work2, double *q, double *rcoef
);
