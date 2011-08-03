/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/tolerance.h"
#include "cskernel/Ble.h"
#include "cskernel/bvabs.h"
#include "cskernel/bpval2.h"
#include "cskernel/Blunk.h"
#include "cskernel/bkcrng.h"
#include "cskernel/blcbp.h"
#include "cskernel/blcpb.h"

#define MAX_ORDER 10
//    BLUCON WILL CONNECT TWO B-REP. CURVES AND GET ONE B-REP. 
// *** INPUT  * 
//     K,N1,T1(N1+K),RCOEF1(IRC1,NCD),IRC1,NCD....1ST B-REP 
//      *** N1 MAY BE 0, IN THIS CASE B-REP 2 COPIED INTO (N1,T1,RCOEF1) 
//     N2,T2(N2+K),RCOEF2(IRC2,NCD),IRC2...2ND B-REP TO CONCATENATE. 
//     JCNI.....C(JCNI) CLASS ABOUT CONTINUITY 
// *** OUTPUT * 
//     N1,T1(N1+K),RCOEF1(.,.).....NEW B-REP CONCATENATED. 
//     RATIO....RATIO OF PARAMETER VALUE REGION (NEW/OLD) 
//                                 OF 2ND B-REP. IN NEW B-REP. 
//     IT2S.....STARTING POINT OF 2ND B-REP. IN NEW B-REP. 
// *** WORK   * 
//     WORK(K,K),PCOEF(K,KK)....
//				WORK ARRAY OF EACH LENGTH, WHERE KK=MAX(NCD,K+2). 
void blucon_(int k, int *n1, double *t1, double *rcoef1, int irc1,
		int ncd, int n2, const double *t2, const double *rcoef2, int irc2, int jcni, 
		double *work, double *pcoef, double *ratio, int *it2s)
{
	double* pcoef1;
    int rcoef1_offset, rcoef2_offset;

    int jcnp1,km1mj, km2mj;
    double delta, d1err, error;
    double rbrk1a[MAX_ORDER], rbrk2a[MAX_ORDER], tworka[MAX_ORDER+MAX_ORDER];
	double *rbrk1, *rbrk2, *twork;

    int i, j, l, m, l1, l2, n3;
    int ib, ii, ki, kk;
    int it, mp1;
    int jcn, n1p1, kpi, nni, n2mk;
    double te,ts,f1,f2, x, x1, x2;
	int modify;

// *********** START OF BLUCON ************************** 
    // Parameter adjustments 
    pcoef1 =pcoef-k-1;
    --t1;
    rcoef1_offset = irc1 + 1;
    rcoef1 -= rcoef1_offset;
    rcoef2_offset = irc2 + 1;
    rcoef2 -= rcoef2_offset;
    --t2;

    if(*n1 <= 0){ // CASE OF NO B-REP IN (N1,T1,RCOEF1).
		*n1 = n2;
	    for (i = 1; i <= n2; ++i) {
			t1[i] = t2[i];
			for (j = 1; j <= ncd; ++j) {
				rcoef1[i+j*irc1] = rcoef2[i+j*irc2];
			}
		}
		for (i = n2 + 1; i <= n2 + k; ++i) t1[i] = t2[i];
		return;
	}

//Prepare local array.
	if(k<=MAX_ORDER){
		rbrk1=rbrk1a; rbrk2=rbrk2a; twork=tworka;
	}else{
		rbrk1=(double *)malloc(sizeof(double)*(k)*4);
		rbrk2=rbrk1+(k); twork=rbrk2+(k);
	}

// *** EVALUATION OF 1ST DERIVATIVE AT CONNECTING POINT BY BLE *** 
	error=bzrzro_();
	jcn = k-2;
    if (jcn >jcni) jcn =jcni;
    x1=t1[*n1+1]; x2=t2[k];
	*ratio = 1.f;
    if(jcn > 0){
		ble_(k,*n1,&t1[1],&rcoef1[rcoef1_offset],irc1,ncd,x1,1,pcoef);
		f1 = bvabs_(ncd, 1, pcoef);
		ble_(k,n2,&t2[1],&rcoef2[rcoef2_offset],irc2,ncd,x2,1,pcoef);
		f2 = bvabs_(ncd, 1, pcoef);
		delta=f1-f2;
		if(delta<0.) delta*=-1.;
		if(delta>error){
			// *** CALCULATION OF ' RATIO ' *** 
			d1err = bzmzro_();
			// *** CHECK OF WHETHER F1(F2) IS 0 OR NOT. IF F1(F2)=0,RATIO=1 ****
			if (f1 > d1err && f2 > d1err) *ratio = f2 / f1;
		}
    }
// ***  CHANGE OF KNOT VECTOR *** 

	if(1.-error<= *ratio && *ratio<=1.+error) modify=0; else modify=1;
		// modify is a flag if knot vector modification of t2 is necessary. 
    jcnp1 = jcn + 1;
    if (jcn <= 0) {
//      *** UPDATE END PART OF 1ST B-REP.  ****** 
	for (i = 1; i <= k; ++i) {
	    twork[i-1] = t1[*n1-k+i];
	    twork[i-1+k] = t1[*n1 + 1];
	}
	blunk_(k, k, &t1[*n1 - k + 1], &rcoef1[*n1 - k + 1 + irc1], 
		irc1, ncd, k, twork, k, work, pcoef);
	for (i = 1; i <= k; ++i) {
	    for (j = 1; j <= ncd; ++j) {
		rcoef1[*n1 - k + i + j * irc1] = pcoef1[i + j *	k];
	    }
	    t1[*n1 + i] = t1[*n1 + 1];
	}
	if (jcn == 0) t1[*n1 + k] = twork[k];
//     *** UPDATE START PART OF 2ND B-REP. *************** 
	for (i = 1; i <= k; ++i) {
	    twork[i-1] = t2[k];
	    twork[i-1+k] = t2[k + i];
	}
	blunk_(k,k,&t2[1],&rcoef2[rcoef2_offset],irc2,ncd,
			k,twork,k,work, pcoef);
//     *** GENERATE T1,N1,RCOEF1 *************** 
	for (i = 1; i <= n2; ++i) {
	    if (t1[*n1 + 1] == t2[k]) {
		t1[*n1 + k - jcnp1 + i] = t2[k + i];
	    } else {
		t1[*n1+k-jcnp1+i] = t1[*n1+1] + t2[k+i] - t2[k];
	    }
	}
	for (i = 1; i <= k; ++i) {
	    for (j = 1; j <= ncd; ++j) {
		rcoef1[*n1 - jcnp1 + i + j * irc1]
			= pcoef1[i + j * k];
	    }
	}
	for (i = 1; i <= n2 - k; ++i) {
	    for (j = 1; j <= ncd; ++j) {
			rcoef1[*n1 + k - jcnp1 + i + j * irc1]
				=rcoef2[k + i + j * irc2];
	    }
	}
	*it2s = *n1 - jcn;
	*n1 = *n1 + n2 - jcnp1;
    } else {
//    **** JCN >= 1  **************************************** 
		ts = t1[*n1 + 1];
		if(modify) te = t1[*n1 + 1] + *ratio * (t2[k + 1] - t2[k]);
		else te = t1[*n1 + 1] + (t2[k + 1] - t2[k]);
		for (i = 1; i <= k; ++i) {
		    twork[i-1] = t2[i];
		    twork[i-1+k] = t2[k + i];
		}
		bkcrng_(k, k, twork, ts, te,twork);
//    **** CHANGE FROM B-REP. TO PP-REP. BY 'BLCBP' 
//              (ONLY CONNECTING SECTION)     **** 
		l = k - jcn;
		for (j = 1; j <= ncd; ++j) {
			blcbp_(k, k, &t1[*n1 - k + 1], &rcoef1[*n1 - k + 1 + j * 
				irc1], work, rbrk1, pcoef, &l1);
		    blcbp_(k, k, twork, &rcoef2[j*irc2 + 1], work, rbrk2,
				&pcoef1[l*k+1], &l2);
		    ib = 1;
//     ***INSERT OF MULTIPLE BREAK POINT SEQUENCE WHICH DEPENDS ON  JCN** 
		    x = rbrk1[1];
		    km1mj = k - 1 - jcn;
		    for (m = 1; m <= km1mj; ++m) {
				mp1 = m + 1;
				for (i = 1; i <= k; ++i)
				    pcoef1[i+ mp1*k] = pcoef1[i+ l*k];
				if (m != km1mj)
				    pcoef1[k + mp1 * k] = pcoef1[k + k];
			}
		    km2mj = k - 2 - jcn;
			if (km2mj != 0) {
			for (i = 1; i <= km2mj; ++i) {
			    m = i + 1;
			    kk = km2mj - i + 1;
				for (ii = 1; ii <= kk; ++ii) {
				pcoef1[k - ii + m * k] =
					bpval2_(rbrk1,pcoef,k,x,ib,k-ii-1);
			    }
			}
			}
//      *** CHANGE FROM PP-REP. TO B-REP. BY BLCPB *** 
			it = 2;
//      *** CREATION OF  NEW KNOT TO INPUT INTO 'BLCPB' *** 
		    n1p1 = *n1 + 1;
		    if (km2mj != 0) {
				for (i = 1; i <= km2mj; ++i) {
					t1[n1p1+i] = t1[n1p1];
				    rbrk1[i+1] = rbrk1[1];
				}
			}
		    rbrk1[l] = rbrk2[1];
		    if (j <= 1) {
				for (i = 1; i <= n2; ++i) {
				    kpi = k + i;
					ki = n1p1 + km2mj + i;
				    if(modify) t1[ki] = x1 + *ratio * (t2[kpi] - x2);
					else t1[ki] = x1 + (t2[kpi] - x2);
				}
			}
		    blcpb_(rbrk1,pcoef,l,k,1,l + 1,irc1,it,
				&t1[*n1 - k + 1], &rcoef1[*n1-k+1+j*irc1], &n3);
		    *it2s = *n1 - jcn;
		    n2mk = n2 - k;
			if (n2mk > 0) {
				for (i = 1; i <= n2mk; ++i) {
				    nni = *n1 - k + n3 + i;
				    rcoef1[nni+j*irc1] = rcoef2[k+i+j*irc2];
				}
		    }
		}
		*n1 = *n1+n2-1-jcn;
    }
	if(k>MAX_ORDER) free(rbrk1);
}
