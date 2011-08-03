/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/b1bfac.h"
#include "cskernel/b1bslv.h"
#include "cskernel/b2vb.h"

// BLGI2D IS SURFACE VERSION OF BLGINT, I.E. 
//   1. RETURNS THE B-COEFFICIENTS IN ROW-WISE 
//   2. CAN BE USED REPEATEDLY BY INPUTTING Q AND IFLAG=1. 
//  BLGI2D  PRODUCES THE B-SPLINE COEFF.S  RCOEF  OF THE SPLINE OF ORDER 
//   K  WITH KNOTS  T(I), I=1,..., N + K , WHICH TAKES ON THE VALUE 
//   P(I,M) AT     TAU(I), I=1,..., N . 
// ******  I N P U T  ****** 
//  IFLAG.....SPECIFIES BLGI2D TO PRODUCE Q FROM TAU(DATA POINT) AND 
//        T(KNOT) -- IFLAG=2 --, OR Q CAN BE UTILIZED BECAUSE 
//        THE PREVIOUS CALL OF BLGI2D PRODUCED THEM AND THIS CALL IS 
//        THE SECOND CALL WITH THE SAME KNOT AND DATA POINT 
//                -- IFLAG=1 --  . 
//  TAU(N)..ARRAY OF LENGTH  N , CONTAINING DATA POINT ABSCISSAE. 
//    A S S U M P T I O N . . .  TAU  IS STRICTLY INCREASING 
//  P(IP,M)...CORRESPONDING ARRAY OF LENGTH  N , CONTAINING DATA POINT 
//       ORDINATES 
//  T(N+K)...KNOT SEQUENCE, OF LENGTH  N+K 
//  N.....NUMBER OF DATA POINTS AND DIMENSION OF SPLINE SPACE 
//  K.....ORDER OF SPLINE 
//  M...NUMBER OF DATA SET(OR SPACE DIMENSION OF POINTS P) 
//  IP....ROW DIMENSION OF THE VARIABLE P 
//  IRC...ROW DIMENSION OF THE VARIABLE RCOEF 
// ******  O U T P U T  ****** 
//  Q.....ARRAY OF SIZE  (2*K-1)*N , CONTAINING THE TRIANGULAR FACTORIZ- 
//       ATION OF THE COEFFICIENT MATRIX OF THE LINEAR SYSTEM FOR THE B- 
//        COEFFICIENTS OF THE SPLINE INTERPOLANT. 
//           THE B-COEFFS FOR THE INTERPOLANT OF AN ADDITIONAL DATA SET 
//        (TAU(I),HTAU(I)), I=1,...,N  WITH THE SAME DATA ABSCISSAE CAN 
//        BE OBTAINED WITHOUT GOING THROUGH ALL THE CALCULATIONS IN THIS 
//       ROUTINE, SIMPLY BY LOADING  HTAU  INTO  RCOEF  AND THEN EXECUT- 
//        ING THE    CALL B1BSLV ( Q, 2*K-1, N, K-1, K-1, RCOEF ) 
//  RCOEF(IRC,N)..THE B-COEFFICIENTS OF THE INTERPOLANT. COEFFICIENTS ARE 
//       STORED IN ROW-WISE OF LENGTH  N. 
//  IFLAG.....AN INTEGER INDICATING SUCCESS (= 1)  OR FAILURE (= 2) 
//       THE LINEAR SYSTEM TO BE SOLVED IS (THEORETICALLY) INVERTIBLE IF 
//        AND ONLY IF 
//              T(I) .LT. TAU(I) .LT. TAU(I+K),    ALL I. 
//        VIOLATION OF THIS CONDITION IS CERTAIN TO LEAD TO  IFLAG = 2 . 
//
// ****** WORK ARRAY ****** 
// WORK(N).....OF LENGTH N 
//
// ******  M E T H O D  ****** 
//    THE I-TH EQUATION OF THE LINEAR SYSTEM  A*RCOEF = B  FOR THE B-CO- 
// EFFS OF THE INTERPOLANT ENFORCES INTERPOLATION AT  TAU(I), I=1,...,N. 
//  HENCE,  B(I) = P(I), ALL I, AND     A  IS A BAND MATRIX WITH  2K-1 
//   BANDS (IF IT IS INVERTIBLE). 
//    THE MATRIX  A  IS GENERATED ROW BY ROW AND STORED, DIAGONAL BY DI- 
//  AGONAL, IN THE  R O W S  OF THE ARRAY  Q, WITH THE MAIN DIAGONAL GO- 
//  ING INTO ROW  K .  SEE COMMENTS IN THE PROGRAM BELOW. 
//     THE BANDED SYSTEM IS THEN SOLVED BY A CALL TO  B1BFAC (WHICH CON- 
//  STRUCTS THE TRIANGULAR FACTORIZATION FOR  A  AND STORES IT AGAIN IN 
//   Q ), FOLLOWED BY A CALL TO  B1BSLV (WHICH THEN OBTAINS THE SOLUTION 
//   RCOEF  BY SUBSTITUTION). 
//    B1BFAC  DOES NO PIVOTING, SINCE THE TOTAL POSITIVITY OF THE MATRIX  
//  A  MAKES THIS UNNECESSARY. 
void blgi2d_(int *iflag, const double *tau, const double *p, 
	const double *t, int k, int n, int m, int ip, 
	int irc, double *work, double *q, double *rcoef){
	int left, lenq;
    double taui;
    int kpkm2, i, j;
    int ilp1mx, jj, km1, np1;

    // Parameter adjustments 
    --work;
    --t;
    --tau;
    p -= ip + 1;;
    rcoef -= irc + 1;
    --q;

    // Function Body 
    km1 = k-1;
    if(*iflag==2){
	    np1 = n + 1;
		kpkm2 = km1 << 1;
	    left = k;
		//                ZERO OUT ALL ENTRIES OF Q 
	    lenq = n*(k+km1);
		for(i=1; i<=lenq; ++i)
			q[i] = 0.f;

		//  *** LOOP OVER I TO CONSTRUCT THE N INTERPOLATON EQUATIONS 
		for(i=1; i<=n; ++i){
			taui = tau[i];
			ilp1mx = i+k;
			if(ilp1mx>np1)
				ilp1mx = np1;
			//** Find LEFT in the closed interval (I,I+K-1) such that 
			//      T(LEFT) .LE. TAU(I) .LT. T(LEFT+1) 
			//   Matrix is singular if this is not possible. 
			if(left<i)
				left=i;
			if(taui<t[left]){
				*iflag = 2;
				return;
			}
			while(taui>=t[left+1]){
				++left;
				if(left>=ilp1mx){
					--left;
					if(taui>t[left+1]){
					    *iflag = 2;
						return;
					}
					break;
				}else
					continue;
			}
			// *** The I-th equation enforces interpolation at TAUI, hence 
			// A(I,J)=B(J,K,T) (TAUI), all J.  Only the K entries with J=  
			// LEFT-K+1,...,LEFT actually might be nonzero. These K numbers  
			// are retutned, in RCOEF  (used for temp.storage here), by the  
			// following. 
				b2vb_(k,n+k,&t[1],taui,left,&work[1]);
			//  We therefore want RCOEF(J)=B(LEFT-K+J)(TAUI) to go into 
			//  A(I,LEFT-K+J), i.e., into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since  
			//  A(I+J,J) is to go into Q(I+K,J), all I,J,if we consider Q 
			//  as a two-dimension array, with 2*K-1 rows (see comments in  
			//  B1BFAC).  In the present program, we treat Q as an equivalen t 
			//  one-dimendional array (because of FORTRAN restrictions on 
			//  DIIMENSION statements). We therefore want RCOEF(J) to go into 
			//  ENTRY 
			//       I-(LEFT+J)+2*K+((LEFT+J)-K-1)*(2*K-1) 
			//             = I-LEFT+1+(LEFT-K)*(2*K-1)+ (2*K-2)*J 
			//  of Q. 
			jj = i-left+1+(left-k)*(k+km1);
			for(j=1; j<=k; ++j){
			    jj += kpkm2;
				q[jj] = work[j];
			}
	    }
	
		//        *** Obtain factoriazation of A, stored again in Q. 
		*iflag=b1bfac_(&q[1], k + km1, n, km1, km1)+1;
	    if(*iflag==2)
			return;
    }

	//***Solve A*RCOEF = P by backward substitution. 
    for(j=1; j<=m; ++j){
		for(i=1; i<=n; ++i)
		    work[i] = p[i+j*ip];
		b1bslv_(&q[1],k+km1,n,km1,km1,&work[1]);
		for(i=1; i<=n; ++i)
			rcoef[j+i*irc]=work[i];
    }
}
