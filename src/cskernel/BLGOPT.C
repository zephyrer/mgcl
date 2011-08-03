/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/b1bfac.h"
#include "cskernel/b1bslv.h"
#include "cskernel/b2vb.h"

//The original codes of this program comes from the FORTRAN code SPLOPT of
//"A Practical Guide to Splines" by Carl de Boor.


#define NEWTON_MAX 10
// blgopt_ COMPUTES THE KNOTS T FOR THE OPTIMAL RECOVERY SCHEME OF ORDER  K 
//  FOR DATA AT  TAU(I), I=1,...,N . 
//
// ******  I N P U T  ****** 
//  TAU.....ARRAY OF LENGTH  N , CONTAINING THE INTERPOLATION POINTS. 
//    A S S U M E D  TO BE NONDECREASING, WITH TAU(I).LT.TAU(I+K),ALL I. 
//  N.....NUMBER OF DATA POINTS . 
//  K.....ORDER OF THE OPTIMAL RECOVERY SCHEME TO BE USED . 
//
// ******  W O R K  A R R A Y  ***** 
//  SCRTCH.....ARRAY OF LENGTH  (N-K)(2K+3) + 5K + 3 . THE VARIOUS 
//        CONTENTS ARE SPECIFIED IN THE TEXT BELOW . 
//
// ******  O U T P U T  ****** 
//  IFLAG.....INTEGER INDICATING SUCCESS (=1) OR FAILURE (=2) . 
//     IF IFLAG = 1, THEN 
//  T.....ARRAY OF LENGTH  N+K  CONTAINING THE OPTIMAL KNOTS READY FOR 
//        USE IN OPTIMAL RECOVERY. SPECIFICALLY,  T(1) = ... = T(K) = 
//        TAU(1)  AND  T(N+1) = ... = T(N+K) = TAU(N) , WHILE THE  N-K 
//        INTERIOR KNOTS  T(K+1), ..., T(N)  ARE CALCULATED AS DESCRIBED 
//        BELOW UNDER  *METHOD* . 
//     IF IFLAG = 2, THEN 
//        K .LT. 3, OR N .LT. K, OR A CERTAIN LINEAR SYSTEM WAS FOUND TO 
//        BE SINGULAR. 
//
// ******  M E T H O D  ****** 
//    THE (INTERIOR) KNOTS  T(K+1), ..., T(N)  ARE DETERMINED BY NEWTONS 
//  METHOD IN SUCH A WAY THAT THE SIGNUM FUNCTION WHICH CHANGES SIGN AT 
//   T(K+1), ..., T(N)  AND NOWHERE ELSE IN (TAU(1),TAU(N)) IS ORTHOGON- 
//  AL TO THE SPLINE SPACE  SPLINE( K , TAU )  ON THAT INTERVAL . 
//     LET  XI(J)  BE THE CURRENT GUESS FOR  T(K+J), J=1,...,N-K. THEN 
//  THE NEXT NEWTON ITERATE IS OF THE FORM 
//              XI(J)  +  (-)**(N-K-J)*X(J)  ,  J=1,...,N-K, 
//  WITH  X  THE SOLUTION OF THE LINEAR SYSTEM 
//                        C*X  =  D  . 
//  HERE,  C(I,J) = B(I)(XI(J)), ALL J, WITH  B(I)  THE I-TH B-SPLINE OF 
//  ORDER  K  FOR THE KNOT SEQUENCE  TAU , ALL I, AND  D  IS THE VECTOR 
//  GIVEN BY D(I) = SUM( -A(J) , J=I,...,N )*(TAU(I+K)-TAU(I))/K, ALL I, 
//  WITH  A(I) = SUM ( (-)**(N-K-J)*B(I,K+1,TAU)(XI(J)) , J=1,...,N-K ) 
//  FOR I=1,...,N-1, AND  A(N) = -.5 . 
//    (SEE CHAPTER  XIII  OF TEXT AND REFERENCES THERE FOR A DERIVATION) 
//    THE FIRST GUESS FOR  T(K+J)  IS  (TAU(J+1)+...+TAU(J+K-1))/(K-1) . 
//     ITERATION TERMINATES IF  MAX(ABS(X(J))) .LT. T O L  , WITH 
//                 T O L  =  T O L R T E *(TAU(N)-TAU(1))/(N-K) , 
//  OR ELSE AFTER  N E W T M X  ITERATIONS , CURRENTLY, 
//                 NEWTMX, TOLRTE / 10, .000001 
//
//     DIMENSION SCRTCH((N-K)*(2*K+3)+5*K+3), T(N+K) 
// URRENT FORTRAN STANDARD MAKES IT IMPOSSIBLE TO SPECIFY THE PRECISE 
// DIMENSIONS OF SCRTCH AND T WITHOUT THE INTRODUCTION OF OTHERWISE 
//  SUPERFLUOUS ADDITIONAL ARGUMENTS . 
void blgopt_(const double *tau, int n, int k, double *scrtch, double *t, int *iflag){
	double tolrte;
    int left, lenw;
    double temp;
    int kpkm1, i, j, l;
    double asign;
    int index, llmin, llmax;
    int id, na, nb, nc, nd, ll, nx;
    double delmax, floatk;
    int leftmk;
    double signst;
    int km1, newton, kp1;
    double del;
    int kpk, nmk, kpn;
    double xij, tol, sum;

	// Parameter adjustments 
    --tau;
    --scrtch;
    --t;

    // Function Body 
    nmk = n-k;
	if(k<=2 || nmk<0){
	    *iflag = 2;
		return;
    }
	
	if(nmk==1){
	    for(i=1; i<=k; ++i){
			t[i] = tau[1];
			t[n+i] = tau[n];
		}
		return;
	}

	tolrte=bzrzro_();
    floatk = (double)(k);
    kpk = k+k;
    kp1 = k+1;
    km1 = k-1;
    kpkm1 = kpk-1;
    kpn = k+n;
    signst = -1.f;
    if(nmk>nmk/2<<1)
		signst = 1.f;

    nx = n + kpk + 1;	//SCRTCH(I) = TAU-EXTENDED(I), I=1,...,N+K+K 
    na = nx + nmk + 1;	//SCRTCH(I+NX) = XI(I),I=0,...,N-K+1 
    nd = na + n;		//SCRTCH(I+NA) = -A(I), I=1,...,N 
    nb = nd + nmk;		//SCRTCH(I+ND) = X(I) OR D(I), I=1,...,N-K 
    nc = nb + kp1;		//SCRTCH(I+NB) = BIATX(I),I=1,...,K+1 
    lenw = kpkm1*nmk;//SCRTCH(I+(J-1)*(2K-1)+NC)=W(I,J)=C(I-K+J,J), I=J-K,...,J+K, J=1,...,N-K.

	//  EXTEND  TAU  TO A KNOT SEQUENCE AND STORE IN SCRTCH. 
    for(j=1; j<=k; ++j){
		scrtch[j] = tau[1];
		scrtch[kpn+j] = tau[n];
    }
    for(j=1; j<=n; ++j)
		scrtch[k+j] = tau[j];

	//  FIRST GUESS FOR  SCRTCH (.+NX)  =  XI . 
    scrtch[nx] = tau[1];
    scrtch[nmk+1+nx] = tau[n];
    for(j=1; j<=nmk; ++j){
		sum = 0.f;
		for(l=1; l<=km1; ++l)
			sum += tau[j+l];
		scrtch[j+nx] = sum/(double)km1;
    }
	//  LAST ENTRY OF  SCRTCH (.+NA)  =  - A  IS ALWAYS ... 
    scrtch[n+na] = .5f;

//  START NEWTON ITERATION. 
// COMPUTE THE 2K-1 BANDS OF THE MATRIX C AND STORE IN SCRTCH(.+NC), 
//  AND COMPUTE THE VECTOR  SCRTCH(.+NA) = -A. 
    tol = tolrte*(tau[n]-tau[1])/(double)nmk;
	for(newton=1; newton<=NEWTON_MAX; newton++){
	    for(i=1; i<=lenw; ++i)
			scrtch[i+nc] = 0.f;
	    for(i=2; i<=n; ++i)
			scrtch[i-1+na] = 0.f;
	    asign = signst;
		left = kp1;
	    for(j=1; j<=nmk; ++j){
			xij = scrtch[j+nx];
			while(xij>=scrtch[left+1]){
				++left;
				if(left>=kpn){
					--left;
					break;
				}
			}
			b2vb_(k, kpn, &scrtch[1], xij, left, &scrtch[nb+1]);
			//THE TAU SEQUENCE IN SCRTCH IS PRECEDED BY  K  ADDITIONAL KNOTS 
			//THEREFORE,  SCRTCH(LL+NB)  NOW CONTAINS  B(LEFT-2K+LL)(XIJ) 
			//WHICH IS DESTINED FOR  C(LEFT-2K+LL,J), AND THEREFORE FOR 
			//    W(LEFT-K-J+LL,J)= SCRTCH(LEFT-K-J+LL + (J-1)*KPKM1 + NC)  
			//SINCE WE STORE THE 2K-1 BANDS OF  C  IN THE 2K-1  R O W S  OF 
			//THE WORK ARRAY W, AND  W  IN TURN IS STORED IN  S C R T C H  , 
			//WITH  W(1,1) = SCRTCH(1 + NC) . 
			//    ALSO, C  BEING OF ORDER  N-K, WE WOULD WANT  1 .LE. 
			//LEFT-2K+LL .LE. N-K  OR 
			//   LLMIN = 2K-LEFT  .LE.  LL  .LE.  N-LEFT+K = LLMAX . 
			leftmk = left-k;
			index = leftmk-j+(j-1)*kpkm1+nc;
			llmin = k-leftmk;
			if(llmin<1)
				llmin = 1;
			llmax = n-leftmk;
			if(llmax>k)
			    llmax = k;
			for(ll=llmin; ll<=llmax; ++ll)
			    scrtch[ll+index] = scrtch[ll+nb];
			b2vb_(kp1,kpn+1,&scrtch[1],xij,left,&scrtch[nb+1]);
			id = leftmk-kp1;
			if(id<0)
				id = 0;
			llmin = 1-leftmk+kp1;
			if(llmin<1)
			    llmin = 1;
			for(ll=llmin; ll<=kp1; ++ll){
				++id;
			    scrtch[id+na] -= asign*scrtch[ll+nb];
			}
			asign = -asign;
	    }
		*iflag=b1bfac_(&scrtch[nc+1],kpkm1,nmk,km1,km1)+1;
	    if(*iflag!=1)
			return;

		// COMPUTE  SCRTCH (.+ND) =  D  FROM  SCRTCH (.+NA) = - A . 
		for(i=n; i>1; i--)
			scrtch[i-1+na] += scrtch[i+na];
	    for(i=1; i<=nmk; ++i)
			scrtch[i+nd] = scrtch[i+na]*(tau[i+k]-tau[i])/floatk;
		// COMPUTE  SCRTCH (.+ND) =  X . 
		b1bslv_(&scrtch[nc + 1], kpkm1, nmk, km1, km1, &scrtch[nd + 1]);
		// COMPUTE  SCRTCH (.+ND) = CHANGE IN  XI . MODIFY, IF NECESSARY, TO 
		//  PREVENT NEW  XI  FROM MOVING MORE THAN 1/3 OF THE WAY TO ITS 
		//  NEIGHBORS. THEN ADD TO  XI  TO OBTAIN NEW  XI  IN SCRTCH(.+NX). 
		delmax = 0.f;
	    asign = signst;
		for(i=1; i<=nmk; ++i){
			del = asign * scrtch[i + nd];
			if(delmax<fabs(del))
				delmax = fabs(del);
			if(del>0.f){
				temp = (scrtch[i+1+nx]-scrtch[i+nx])/3.f;
			    if(del>temp)
					del = temp;
			}else{
			    temp = (scrtch[i-1+nx]-scrtch[i+nx])/3.f;
				if(del<temp)
					del = temp;
			}
			asign = -asign;
			scrtch[i+nx] += del;
		}		
		if(delmax<tol)//IN CASE CHANGE IN XI WAS SMALL ENOUGH OR TOO MANY STEPS WERE TAKEN.
			break;
	}

    for(i=1; i<=nmk; ++i)
		t[k+i] = scrtch[i+nx];
    for(i=1; i<=k; ++i){
		t[i] = tau[1];
		t[n+i] = tau[n];
    }
}
