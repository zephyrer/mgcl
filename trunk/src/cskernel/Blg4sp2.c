/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/b1bfac.h"
#include "cskernel/b1bslv.h"
#include "cskernel/b2vb.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

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
){
    int nend, left, lenq;
    double taui;
    int kpkm2, i, j;
    int iq, lftmax;
    double dt1, dt2;
	int kpkm1, km1, np1,npk;
    double dt12;
    double fkm1;
// ****************** START OF BLG4SP ***********************************
    // Parameter adjustments 
    --tau;
    --work2;
    --work1;
    val -= iv+1;
    --t;
    rcoef -= irc+1;
    --q;

    km1 = k-1;
    kpkm1 = k+km1;
	if(*iflag == 2){

    if(ibcbeg<=0 || ibcbeg>4 || ibcend<=0|| ibcend>4)
		goto L7000;

    if(n<k)
		goto L7000;

//   ********* THE B-REP ORDER K IS ALWAYS 4 ******** 
    fkm1 = (double) km1;
    kpkm2 = kpkm1-1;
    np1 = n+1;
    npk = n+k;
    if(t[k+1] <= t[k] || t[np1] <= t[n])
		goto L7000;

//    CLEAR ALL ENTRIES OF Q(.) 
    lenq = n * kpkm1;
    for (i = 1; i <= lenq; ++i)
		q[i] = 0.f;
 
//    FIRST POINT MAKES FIRST B-COEFFICIENT. 
    i=1;
    if(tau[i] != t[k])
		goto L7000;

	iq=k;
    left=k;
    q[iq] = 1.f;
    work2[i] = 1.f;
// =====SET B.C. OF IBCBEG.===== 
    ++iq;
	dt1 = t[k+1]-t[k];
	if(ibcbeg==1 || ibcbeg==4){
	//  1ST DERIV. SPECIFIED IN VAL(2,.) 
	    ++i;
		work2[i] =dt1/(fkm1*2.f);
	    q[iq] = -.5f;
		q[iq+kpkm2] = .5f;
	    ++iq;
	}
	if(ibcbeg==2 || ibcbeg==4){
	//  2ND DERIV. SPECIFIED. 
	    ++i;
	    dt2 = t[k+2]-t[k];
		dt12 = (dt1+dt2)*2.f;
	    q[iq] = dt2/dt12;
		q[iq+kpkm2] = -.5f;
	    q[iq+(kpkm2<<1)] = dt1/dt12;
		work2[i] = dt1*dt1*dt2/(fkm1*(double)(k-2)*dt12);
	}
	if(ibcend==1 || ibcend== 2)
		nend = n-2;
	else if(ibcend==3)
	    nend = n-1;
	else
		nend=n-3;

// =====LOOP OVER I TO CONSTRUCT THE INTERPOLANT EQUATIONS.===== 
	for(++i;i<=nend;i++){
	    taui = tau[i];
		//       CLARIFY WETHER MULTIPLE DATA POINT OR NOT 
		if(i==nend || taui!=tau[i+1]){
		// CASE OF NO DATA POINT MULTIPLICITY, LOCATE WHERE TAUI IN T(.) 
		    lftmax = i + k;
		    if(lftmax>np1)
				lftmax=np1;
		    if(left<i)
				left=i;
		    if(taui<t[left])
				goto L7000;
	
			while(taui>=t[left+1]){
				++left;
				if(left>=lftmax)
					goto L7000;
			}
			//NOW T(LEFT) <= TAU(I) < T(LEFT+1), GET B-SPLINE BASIS FUNCTION AND STORE IN Q(.) 
			b2vb_(k, npk, &t[1], taui, left, &work1[1]);
			iq = i-left+1+(left-k)*kpkm1;
		    for(j=1; j<=k; ++j){
				iq += kpkm2;
				q[iq] = work1[j];
		    }
			work2[i] = 1.f;
		}else{
		//CASE OF MULTIPLE DATA POINT 1, 1ST DERIV. SPECIFIED. 
		    if(taui!=t[i+2] || taui!=t[i+3])
				goto L7000;

			// ..FOR 1ST DERIV. DATA 
		    iq = i*kpkm1-km1;
		    q[iq] = -.5f;
		    q[iq+kpkm2] = .5f;
		    work2[i] = (t[i+4]-t[i+1])/(fkm1*2.f);

			//      CHECK 2ND MULTIPLICITY. 
		    if(i+2>nend)
			    continue;
		    if(taui<tau[i+2])
			    continue;
	
			if(taui!=t[i+4])
				goto L7000;

			//      CASE OF MULTIPLE DATA POINT 2. 
			//                1ST DERIV. DISCONTINUITY. 
			//         ..FOR POSITIONAL DATA(POSITIONAL DATA MAKE B-COEF) 
			++i;
		    iq = i*kpkm1-km1;
		    q[iq] = 1.f;
		    work2[i] = 1.f;
			//         ..FOR 1ST DERIV DATA 
		    ++i;
		    ++iq;
		    q[iq] = -.5f;
		    q[iq+kpkm2] = .5f;
		    work2[i] = (t[i+3]-t[i])/(fkm1*2.f);
		}
	}

// =====   SET B.C. OF IBCEND ===== 
    if (tau[n] != t[np1])
		goto L7000;
    iq = i*kpkm1-km1;	
	if(ibcend==4){
	    iq += kpkm2;
	}
    dt1 = t[np1]-t[n];
	if(ibcend==2 || ibcend==4){
		dt2 = t[np1]-t[n-1];
	    dt12 = (dt1+dt2)*2.f;
		q[iq-kpkm2] = dt1/dt12;
	    q[iq] = -.5f;
		q[iq+kpkm2] = dt2/dt12;
	    work2[i] = dt1*dt1*dt2/(fkm1*(double)(k-2)*dt12);
		if(ibcend==4){
		    ++iq; ++i;
		}
	}
	if(ibcend==1 || ibcend==4){
	    q[iq] = -.5f;
		q[iq + kpkm2] = .5f;
	    work2[i] = dt1/(fkm1*2.f);
	}

    q[n*kpkm1-km1] = 1.f;
    work2[n] = 1.f;
// ===== OBTAIN FACTORIZATION OF THE MATRIX, STORED AGAIN IN Q.===== 
    *iflag=b1bfac_(&q[1], kpkm1, n, km1, km1)+1;
    if(*iflag == 2)
		return;

	}

// =====  SOLVE A*RCOEF = VAL. ===== 
    for(j=1; j<=m; ++j){
		for(i=1; i<=n; ++i){
			work1[i] = work2[i]*val[i+j*iv];
		}
		b1bslv_(&q[1], kpkm1, n, km1, km1, &work1[1]);
		for(i=1; i<=n; ++i)
			rcoef[j+i*irc] = work1[i];
	}
    return;

// ==== ERROR RETURN ===== 
L7000:
    *iflag = 2;
}
