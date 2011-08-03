/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkktdp.h"

// BLUMOR         BLUMOR FOR BLUMOV 
// BLUMOR GENERATES RATIO'S OF B-COEFFICIENT TRANSLATION. 
// *** INPUT * 
//       KMOVI.....INDICATES WHAT KIND OF MOVE IS BEING PERFORMED. 
//               =1 : TWO POINTS FIXED AND CENTER OF MOVE BETWEEN THEM. 
//               =2 : ONE POINT,T(I),FIXED. THE OTHER POINT IS FREE 
//               =3 : ONE POINT,T(J+1),FIXED. THE OTHER POINT IS FREE 
//       I,J....GIVE ID OF KNOT VECTOR T(.) THAT SHOULD BE FIXED. 
//                 T(I)<= TAU <=T(J+1)   K<=I<=N, K<=J<=N. 
//       TAU......IS THE CENTER OF TRANSLATION, PARAMETER VALUE OF THE 
//                B-REP. 
//       K,N,T....DESCRIBE KNOT VECTOR OF THE B-REP. 
//                K : ORDER,      N : B-REP DIMENSION 
//                T(N+K) : KNOT VECTOR 
// *** OUTPUT * 
//       I,J.....UPDATED I AND J. 
//       ISTRT,IEND......START AND END ID OF B-COEF RATIO(.,.) THAT 
//                SHOULD BE MODIFIED. 
//       RATIO(M)...RATIO(M) CONTAINS RATIO OF TRANSLATION 
//                  RATIO(M) IS FOR RCOEF(M,.) 
//                    ISTRT<= M <= IEND 
// ***WORK* 
//       RATIO(N)...WORK AREA OF LENGTH N. 
void blumor_(int kmovi, int *i, int *j, 
	double tau, int k, int n, const double *t, int *istrt,
	int *iend, double *ratio)
{
    // Local variables 
    int kmov;
    double tlen1, tlen2;
    int m, i1, i2;
    double t1, t2, taumt1, t2mtau;
    int np1;
    double t1l, t2l;

   // Parameter adjustments 
    --ratio;
    --t;

    // Function Body 
    kmov = kmovi;
    np1 = n+1;

    t1 = t[*i];
    taumt1 = tau-t1;
    t1l = (t[*i + 1] - t[*i]) / 3.f;
    if(taumt1<t1l){
		--(*i);
		while(*i>k && t[*i-1] == t[*i])
			--(*i);
		if(*i<=k && t[k]==t1){
		    kmov = 3;
			*i = k;
		}
    }

    t2 = t[*j+1];
    t2mtau = t2-tau;
    t2l = (t[*j+1]-t[*j])/3.f;
    if(t2mtau<t2l){
		++(*j);
		while(*j<n && t[*j]==t[*j+1])
			++(*j);
		if(*j>=n && t[np1]==t2){
			kmov = 2;
		    *j = n;
		}
    }

    if(kmov==1 || kmov==2){
		*istrt = *i-k+2;
    }else{
		*istrt = 1;
    }
    if(kmov==1 || kmov==3){
		*iend = *j-1;
    }else{
		*iend = n;
    }

    bkktdp_(n, k, &t[1], &ratio[1]);
    i1 = *istrt-1;
    if(i1<1)
		i1 = 1;
    i2 = *iend + 1;
    if(i2>n)
		i2 = n;
    t1 = ratio[i1];
    t2 = ratio[i2];
    tlen1 = tau-t1;
    tlen2 = t2-tau;
    if(kmov==1){
		for(m=*istrt; m<=*iend; ++m){
			if(ratio[m]<=tau){
				ratio[m] = (ratio[m]-t1)/tlen1;
		    } else {
				ratio[m] = (t2-ratio[m])/tlen2;
		    }
		}
    }else if(kmov==2){
		for(m=*istrt; m<=*iend; ++m){
		    ratio[m] = (ratio[m]-t1)/tlen1;
		}
    }else{
		for(m=*istrt; m<=*iend; ++m)
			ratio[m] = (t2-ratio[m])/tlen2;
    }
}
