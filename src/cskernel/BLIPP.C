/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bler.h"
#include "cskernel/blip1.h"
#include "cskernel/bli1sp.h"
#include "cskernel/bk1fli.h"
#include "cskernel/blipp3.h"

// BLIPP COMPUTES INTERSECTION PARAMETER VALUE X(I) OF 1-D B-REP, S.T. 
//  F = G(X(I)) FOR G(.):B-REP FUNCTION, AND X(I) > TS FOR 1<=I<=MX. 
// *** INPUT * 
//     K,N,T(N+K),RCOEF(N)......BREP OF ONE SPACE DIMENSION. 
//     F......  B-REP FUNCTION VALUE TO GET THE PARAM VALUES. 
//     ERROR.....ERROR ESTIMATE ALLOWED TO GET THE INTERSECTION POINTS. 
//     TS.......INDICATES PARAMETER VALUE AT WHICH TO START THE 
//              COMPUTATION. 
//     MX......PROVIDES LENGTH OF THE VARIABLE X(.), BLIPP STOPS 
//             COMPUTATION WHEN NUM OF INTERSECTION POINTS REACHES MX. 
// *** OUTPUT * 
//     NX,X(NX)..... PARAMETER VALUES OBTAINED BY BLIPP. NX<=MX. 
//                   NX=0 WHEN NO SOLUTION. 
//     IEND.....GIVES THE INF WETHER BLIPP COMPUTES TO THE END OF 
//              PARAMETER. 
//               =0: NOT TO THE END BECAUSE MX IS TOO SMALL. 
//               =1: COMPUTE TO THE END. 
// *** WORK * 
//     WORK(.)       WORK ARRAY OF LENGTH 4*K*K+3*K 
void blipp_(int k, int n, const double *t, const double *rcoef,
	double f, double error, double ts, int mx,
	double *work, int *nx, double *x, int *iend)
{
    double d;
    int isgn;
    double tnew;
    int iend3;
    int isold, nspan, isnew;
    double d1, d2;
    double fm;
    int ki, it;
    double tk;
    int km2, itnext;
    double errort;
    int np1, iw2, iw3, ith, itl, kip1, kby2;
    double tnp1;

// **************************** START OF BLIPP************************* 
    // Parameter adjustments 
    --rcoef;
    --t;
    --x;

    // Function Body 
    *nx = 0;
    *iend = 0;
    if (mx <= 0)
		return;

    np1 = n+1;
    tnp1 = t[np1];
	if(ts>=tnp1){
		*iend = 1;
	    return;
	}

    km2 = k-2;
    iw2 = k*k+k;
    iw3 = k*3*k+k*3;
    kby2 = k/2;
    tk = t[k];

	errort=(tnp1-tk)*bzrzro_();
	errort*=5.;					// Updated on 1999/7/29. by mizuno.

    if(ts>=tk){
		ki=bk1fli_(np1, &t[1], ts);
		ki=ki-k+1;
    }else{
	// CHECK OF THE STARTING POINT. 
		ki = 1;
		fm = bler_(k,n,&t[1],&rcoef[1],tk,0);
		if(fabs(fm-f)<=error){
		    *nx=1;
			x[1]=tk;
		}
    }

// ++++LOOP OVER KI UNTIL KI => N . 
	while(1){
		d = rcoef[ki]-f;
		if(d<0.)
			isold=-1;
	    else if(d>0.)
			isold = 1;
	    else
			isold = 0;

		while(1){
			if(ki>=n)
				break;
		    kip1=ki+1;
			d=rcoef[kip1]-f;
			if(d<0.)
				isnew = -1;
		    else if(d>0.)
				isnew = 1;
			else
				isnew = 0;
			if(isold==0 || isold*isnew<0)
				break;
			isold=isnew;
			ki=kip1;
		}
		if(ki>=n)
			break;

		//  NOW  F IS IN BETWEEN RCOEF(KI),RCOEF(KI+1), 
		//  CHECK THE INTSCTN NUM. OF THE NEXT K-2 SPANS. 
	    nspan = n-kip1;
		if(nspan>km2)
			nspan=km2;

		if(nspan==0 || blip1_(&rcoef[1],kip1,nspan,f)!=0){
		// CASE OF COMPLICATED SOLUTION, OBTAIN THE SOLUTIONS BY BLIPP3. 
		    blipp3_(k,n,&t[1],&rcoef[1],&ki,f,error,errort,mx,ts,work, 
				work+iw2, work+iw3, nx, &x[1], &iend3);
		    if(iend3==0)
				return;
			continue;
		}

		//     NOW NUM OF SOLUTIONS IS AT MOST ONE IN T(KI+1)<=  <T(KI+K) . 
		itl = kip1;
	    if(itl<k)
			itl=k;
		ith = ki+k;
	    if(ith>np1)
			ith=np1;
		ki = kip1+nspan;

	    d1 = f-bler_(k,n,&t[1],&rcoef[1],t[itl],0);
		if(fabs(d1)>error){
		    d2 = f-bler_(k,n,&t[1],&rcoef[1],t[ith],0);
			if(fabs(d2)>error){
				if(d1*d2>0.f)
					continue;
				// NOW THE ONLY ONE SOLUTION EXIST. 
			    it = itl+kby2;
				if(it>ith)
					it=ith;
				d2 = f-bler_(k,n,&t[1],&rcoef[1],t[it],0);
				isgn = 1;
			    if(d1*d2<0.f)
					isgn=-1;
				while(1){
					itnext = it+isgn;
					d1 = f-bler_(k,n,&t[1],&rcoef[1],t[itnext],0);
					if(d1*d2<= 0.f)
						break;
					d2 = d1;
					it = itnext;
				}
				if(it>itnext)
					it=itnext;
				tnew = bli1sp_(k,n,&t[1],&rcoef[1],it,f,error,work);
			}else
			    tnew = t[ith];
		}else
			tnew = t[itl];
		if(tnew<=ts)
			continue;

		if(*nx==0 || (tnew-errort)>x[*nx])
			++(*nx);
		if(*nx>=mx)
			return;
		x[*nx]=tnew;
	}

// CASE OF THE LAST B-COEFFICIENT. 
    fm = bler_(k,n,&t[1],&rcoef[1],tnp1,0);
	if(fabs(fm-f)<=error){
		if(*nx==0 || (tnp1-errort)>x[*nx])
		    ++(*nx);
	    if(*nx>=mx)
			return;
	    x[*nx] = tnp1;
	}
    *iend = 1;
}
