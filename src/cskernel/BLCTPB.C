/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bpval.h"
#include "cskernel/blcpb2.h"

// BLCTPB COMOUTES B-COEFFICIENTS OF THE B-REP OF PP-REP(RBRK,PCOEF,L). 
// ***INPUT* 
//    RBRK(L+1)........BREAK POINT SEQUENCE OF PP-REP. 
//    PCOEF(K,IPC,NCD).PP-COEFFICIENTS , I.E. 
//             PCOEF(J,I,,)=D**(J-1)(F(RBRK(I))   1<=I<=L. 
//    L................INDICATES NUM OF INTERVAL OF PP-REP. 
//    K................ORDER OF PP-REP(B-REP) 
//    NCD..............SPACE DIMENSION OF THE PP-REP 
//    N,T(N+K).........KNOT VECTOR OF THE BREP. N IS B-REP DIMENSION. 
//    IPC..............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCOEF. 
//    IPCW.............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCWORK. IPCW MUST BE GRETAER OR EQUAL TO (N-K+2). 
//    IRC..............ROW DIMENSION OF THE VARIABLE RCOEF 
// ***OUTPUT*** 
//    RCOEF(IRC,NCD)...B-COEFFICIENTS 
// ***WORK*** 
//    PCWORK(K,IPCW,NCD)........IS WORK ARRAY. IPCW>=(N-K+2). 
void blctpb_(
	const double *rbrk, const double *pcoef, int l,int k,
	int ncd, int n, const double *t, int ipc, int ipcw, int irc,
	double *pcwork, double *rcoef
){
    int pcwork_offset;
    int mold, mmim, mult, nmkp1, nmkp2, i, j, m;
    int im, im1, km1, mm1, lp1, mom1;

    // Parameter adjustments 
    --rbrk;
    --t;
    pcoef -= k*(ipc+1)+1;
    pcwork_offset = k*(ipcw+1)+1;
    pcwork -= pcwork_offset;
    rcoef -= irc+1;

    // Function Body 
    for(j=1; j<=ncd; ++j){
		for(i=1; i<=k; ++i)
			pcwork[i+(j*ipcw+1)*k] = pcoef[i+(j*ipc+1)*k];
    }
    km1 = k-1;
    nmkp1 = n-k+1;
    mold = 1;
    mult = 0;
    for(m=2; m<=nmkp1; ++m){
		mm1=m-1;
		if(t[k+mm1] == t[k+m-2]){
			++mult;
			for(im=1; im<=mult; ++im){
				mmim = m-im;
				for(j=1; j<=ncd; ++j){
				    for(i=1; i<=k; ++i)
						pcwork[i+(mmim+1+j*ipcw)*k] = pcwork[i+(mmim+j*ipcw)*k];
				}
			}
			mom1 = mold-1;
		    if(mult==1){
				for(j=1; j<=ncd; ++j)
				    pcwork[k+(mm1+j*ipcw)*k]=pcoef[k+(mom1+j*ipc)*k];
			}
			for(j=1; j<=ncd; ++j)
				pcwork[k-mult+(m-mult+j*ipcw)*k]
					=bpval_(&rbrk[mom1],&pcoef[(mom1+j*ipc)*k+1],1,k,rbrk[mold],km1-mult);
		}else{
			mult = 0;
		    ++mold;
		    for(j=1; j<=ncd; ++j){
				for(i=1; i<=k; ++i)
					pcwork[i+(m+j*ipcw)*k] = pcoef[i+(mold+j*ipc)*k];
			}
		}
    }
    lp1 = l+1;
    nmkp2 = n-k+2;
    for(i=1; i<=k; ++i){
		im1=i-1;
		for(j=1; j<=ncd; ++j)
			pcwork[i+(nmkp2+j*ipcw)*k]
				=bpval_(&rbrk[1],&pcoef[(j*ipc+1)*k+1],l,k,rbrk[lp1],im1);
    }
	// GENERATE B-COEFFICIENTS 
    blcpb2_(&pcwork[pcwork_offset],nmkp1,k,n,&t[1],ncd,ipcw,irc,&rcoef[irc+1]);
}
