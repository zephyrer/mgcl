/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>

//The original codes of this program comes from the FORTRAN code BSPLPP of
//"A Practical Guide to Splines" by Carl de Boor.

// CONVERTS THE B-REPRESENTATION  K,N,T,BCOEF OF SOME SPLINE INTO ITS 
//  PP-REPRESENTATION  tau, COEF, L, K . 
// ******  I N P U T  ****** 
//  ORDER.....ORDER OF THE SPLINE 
//  N.....LENGTH OF  BCOEF  AND  DIMENSION OF SPLINE SPACE  SPLINE(K,T) 
//  T.....KNOT SEQUENCE, OF LENGTH  N+K 
//  BCOEF(irc,ncd).....B-SPLINE COEFFICIENT SEQUENCE, OF LENGTH  N,
//				and space dimension ncd
//The row dimension of BCOEF is irc.
//ipc.....indicates the structure of coef as described in coef.
// ******  W O R K   A R E A  ****** 
//  SCRTCH......OF SIZE  (K,K,ncd) , NEEDED TO CONTAIN BCOEFFS OF A PIECE OF 
//                THE SPLINE AND ITS  K-1  DERIVATIVES 
// ******  O U T P U T  ****** 
//  TAU.....BREAKPOINT SEQUENCE, OF LENGTH  L+1, CONTAINS (IN INCREAS- 
//        ING ORDER) THE DISTINCT POINTS IN THE SEQUENCE  T(K),...,T(N+1) 
//  COEF(K,ipc,ncd).....ARRAY OF SIZE (K,ipc,ncd), WITH  COEF(I,J, m) =
//		(I-1)ST DERIVATIVE of m-th space dimension OF SPLINE AT tau(J)
//		FROM THE RIGHT                                                  
//  L...NUMBER OF POLYNOMIAL PIECES WHICH MAKE UP THE SPLINE IN THE
//		INTERVAL  (T(K), T(N+1)). L's maximum is N+1-K.                 
// ******  M E T H O D  ****** 
//   FOR EACH BREAKPOINT INTERVAL, THE  K  RELEVANT B-COEFFS OF THE 
//   SPLINE ARE FOUND AND THEN DIFFERENCED REPEATEDLY TO GET THE B-COEFFS 
//   OF ALL THE DERIVATIVES OF THE SPLINE ON THAT INTERVAL. THE SPLINE AND 
//   ITS FIRST  K-1  DERIVATIVES ARE THEN EVALUATED AT THE LEFT END POINT 
//   OF THAT INTERVAL. 
void blcbpn_(int order, int n, const double *t,const double *bcoef, int irc, int ncd, int ipc,
	double *scrtch,	double *tau, double *coef, int *l){
    // Local variables 
    double diff, fkmj;
    int left, leftp1;
    double term;
    int i, j, k, m, js, jsi,jsi1, jsimk, im1;
    int lsofar, lsoc, jp1, kmj, kmjs, kbyk, kipc, i1,i2;
    double sum, saved;
    double biatx_buf[20], deltal_buf[20], deltar_buf[20];
	double* deltal=deltal_buf; double* deltar=deltar_buf;
	double* biatx=biatx_buf;

	k=order;
	if(k>20){
		// type cast (void* -> double*)
		biatx=(double*)malloc(sizeof(double)*(k));
		deltal=(double*)malloc(sizeof(double)*(k));
		deltar=(double*)malloc(sizeof(double)*(k));
	};

    // Parameter adjustments 
    coef -= 1+k+k*ipc;
    scrtch -= 1+k+k*k;
    bcoef-=1+irc;
    --t;
    --tau;

    // Function Body 
	kbyk=k*k;
	kipc=k*ipc;
    lsofar = 0;
    tau[1] = t[k];
    for (left = k; left <= n; ++left){
//                            FIND THE NEXT NONTRIVIAL KNOT INTERVAL.

	leftp1=left + 1;
	if (t[leftp1] == t[left]) continue;
	++lsofar;
	lsoc=lsofar*k;

	tau[lsofar+1] = t[leftp1];
	if(k<=1){
		i1=1+lsoc;
	    for(m=1;m<=ncd;m++) coef[i1+m*kipc] = bcoef[left+m*irc];
		continue;
	}
//     STORE THE K B-SPLINE COEFF.S RELEVANT TO CURRENT KNOT INTERVAL
//                             IN  SCRTCH(.,1) . 
	for(i=1; i<=k; ++i){
		i1=i+k; i2=left-k+i;
	    for(m=1;m<=ncd;m++)
			scrtch[i1+m*kbyk] = bcoef[i2+m*irc];
	}

//     FOR J=1,...,K-1, COMPUTE THE  K-J  B-SPLINE COEFF.S RELEVANT TO
//     CURRENT KNOT INTERVAL FOR THE J-TH DERIVATIVE BY DIFFERENCING 
//     THOSE FOR THE (J-1)ST DERIVATIVE, AND STORE IN SCRTCH(.,J+1). 
	for(jp1=2; jp1<=k; ++jp1){
	    j = jp1-1;
	    kmj = k-j;
	    fkmj = (double)kmj;
		js=j*k;
	    for(i=1; i<=kmj; ++i){
			jsi=i+js; jsi1=i+jp1*k;
			diff = t[left+i] - t[left+i-kmj];
			if(diff>0.f){
				for(m=1;m<=ncd;m++){
					jsimk=jsi+m*kbyk;
				    scrtch[jsi1+m*kbyk]=(scrtch[jsimk+1]-scrtch[jsimk])/diff*fkmj;
				}
			}
	    }
	}

//     FOR  J = 0, ..., K-1, FIND THE VALUES AT  T(LEFT)  OF THE  J+1 
//     B-SPLINES OF ORDER  J+1  WHOSE SUPPORT CONTAINS THE CURRENT 
//     KNOT INTERVAL FROM THOSE OF ORDER  J  (IN  BIATX ), THEN COMBINE
//     WITH THE B-SPLINE COEFF.S (IN SCRTCH(.,K-J) ) FOUND EARLIER 
//     TO COMPUTE THE (K-J-1)ST DERIVATIVE AT  T(LEFT)  OF THE GIVEN 
//     SPLINE. 

	biatx[0]=1.f;
	for(m=1;m<=ncd;m++)
		coef[k+lsoc+m*kipc] = scrtch[1+kbyk+m*kbyk];
	for(jp1=2; jp1<=k; ++jp1){
	    j = jp1-1;
	    deltar[j-1] = t[left+j]-t[left];
	    deltal[j-1] = t[left]-t[leftp1-j];
	    saved = 0.f;
	    for(i=1; i<=j; ++i){
			im1=i-1;
			term = biatx[im1]/(deltar[im1] + deltal[j-i]);
			biatx[im1] = saved+deltar[im1]*term;
			saved = deltal[j-i]*term;
	    }
	    biatx[j] = saved;
	    kmj = k-j; kmjs=kmj*k;
		for(m=1;m<=ncd;m++){
		    sum = 0.f;
		    for(i=1; i<=jp1; ++i){
				sum = biatx[i-1]*scrtch[i+kmjs+m*kbyk]+sum;
		    }
			coef[kmj+lsoc+m*kipc] = sum;
		}
	}

    }

    *l = lsofar;
	if(k>20){ free(biatx); free(deltal); free(deltar);}
}
