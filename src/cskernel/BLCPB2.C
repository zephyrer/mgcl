/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
// BLCPB2 COMPUTES B-COEFFICIENTS OF B-REP, FROM P-REP.(RBRK,PCOEF,L) AND KNOT VECTOR T. 
// ***INPUT* 
//    PCOEF(K,IPC,NCD).PP-COEFFICIENTS , I.E. 
//             PCOEF(J,I,,)=D**(J-1)(F(RBRK(I))   1<=I<=L. 
//    L................INDICATES NUM OF INTERVAL OF PP-REP. 
//    K................ORDER OF PP-REP(B-REP) 
//    N................B-REP DIMENSION 
//    T(N+K)...........KNOT VECTOR OF B-REP 
//    NCD..............SPACE DIMENSION OF THE PP-REP 
//    IPC..............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCOEF, MUST BE .GE.L+1 . 
//    IRC..............ROW DIMENSION OF THE VARIABLE RCOEF 
// ***OUTPUT*** 
//    RCOEF(IRC,NCD)...B-COEFFICIENTS 
void blcpb2_(const double *pcoef, int l, int k, int n,
			const double *t, int ncd, int ipc, int irc, double *rcoef){
    int km2mj, i, j,km2mi2, i1, i2, i3, i4;
    int ip, kd2, km1, km2, lp1, np1, kd2p;
    double r1;

    double x_buf[18], tx_buf[18];
	double* x=x_buf; double* tx=tx_buf;
	if(k>20){
		// type cast (void* -> double*)
		x=(double*)malloc(sizeof(double)*(k-2));
		tx=(double*)malloc(sizeof(double)*(k-2));
	}

    // Parameter adjustments 
    pcoef -= k*(ipc+1)+1;
    rcoef -= irc+1;
    --t;

    // Function Body 
    lp1 = l+1;
    np1 = n+1;
    km1 = k-1;
    km2 = k-2;
    kd2 = k/2;
    kd2p = (k + 1) / 2;
// GENERATE B-COEFFICIENTS (LOOP OVER I UNTIL N) 
    for(i=1; i<=n; ++i){

	ip=i-kd2p+1;
	if(ip<1)
	    ip = 1;
	if(ip>lp1)
	    ip = lp1;
//        IP IS THE INDEX OF PCOEF(.,IP,.) 
	i3=i+kd2;
	if(i3<k)
		i3=k;
	if(i3>np1)
	    i3=np1;
//     GENERATE M.T.I,J    J=1...K-2 
	i1 = 1;
	for(i2=1; i2<=km1; ++i2){
	    i4 = i+i2;
	    if(i4==i3)
			continue;
	    x[i1-1] = t[i4]-t[i3];
	    ++i1;
	}
	for(i2 = 1; i2 <= km2; ++i2)
	    tx[i2-1] = 1.f;
	for(j = 1; j <= km2; ++j) {
	    km2mj = km2-j;
	    tx[km2-1] = x[km2mj]*tx[km2-1];
	    for(i2=1; i2<=km2mj; ++i2){
			km2mi2 = km2-i2;
			tx[km2mi2-1]=tx[km2mi2]+x[km2mj-i2]*tx[km2mi2-1];
	    }
	}
//        NOW M.T.I,J  IN TX(J). 
//  GET BCOEFFICIENTS IN RCOEF(I,.). 
	for(i1=1; i1<=ncd; ++i1){
	    rcoef[i+i1*irc] = 0.f;
	    r1 = 2.f;
	    for(i2=1; i2<=km2; ++i2){
			rcoef[i+i1*irc]=(rcoef[i+i1*irc]+tx[km1-i2-1]*pcoef[k-i2+(ip+i1*ipc)*k])/r1;
			r1 += 1.f;
	    }
	    rcoef[i+i1*irc] += pcoef[(ip+i1*ipc)*k+1];
	}

    }
	if(k>20) {free(x); free(tx);}
}
