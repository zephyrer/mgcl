/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/bkdtpg.h"
#include "cskernel/bkdnp.h"
#include "cskernel/Ble.h"
#include "cskernel/blg4sq.h"

//  BLG4SC GENERATES CIRCULAR B-SPLINE, THAT IS, STARTING AND ENDING 
//  POINTS COINCIDE AND THEIR TANGENT'S ARE EQUAL. 
// *** INPUT * 
//         NV,VAL(IV,NCD),IV,NCD.....ORIGINAL DATA OF SPACE DIMENSION 
//                NCD. NV IS NUM OF DATA AS VAL(I,.) 1<=I<=NV. 
//         IRC....ROW-DIMENSION OF RCOEF AS RCOEF(IRC,NCD). 
// *** OUTPUT * 
//         N,T(N+4),RCOEF(IRC,NCD).....B-SPLINE OF ORDER 4 OBTAINED. 
//                THE SPACE DIMENSION IS NCD AND THE B-REP DIMENSION 
//                IS N=NV+2. 
//         IFLAG  : =1 NORMAL END 
//                  <>1 ABNORMAL 
// *** WORK * 
//         TAU(N),WORK(N,9) : WORK AREA FOR SUBROUTINE BLG4SQ 
void blg4sc_(int nv, double *val, int iv, 
	int ncd, int irc, double *tau, double *work, 
	int *n, double *t, double *rcoef, int *iflag
){
    // System generated locals 
    int val_offset, rcoef_offset, nvm1;

    // Local variables 
    double taul[5],tangn[3];
    int i, j, nm1;

// *** 1. GENARATE DATA POINTS INCLUDING END POINTS' TANGENT IN TAU. 
    // Parameter adjustments 
    val_offset = iv+1;
    val -= val_offset;
    rcoef_offset = irc+1;
    rcoef -= rcoef_offset;
    --tau;
    --t;

    // Function Body 
    bkdtpg_(&val[val_offset], nv, ncd, iv, &tau[2]);
    bkdnp_(&nv, &tau[2], &val[val_offset], iv, ncd, 1, bkmax_());
    tau[1] = tau[2];
    *n = nv + 2;
    nm1 = *n - 1;
    tau[*n] = tau[nm1];
// *** 2. COMPUTE TEMPORARY B-REP TO GET TANGENT, THAT IS, STARTING 
//        AND ENDING POINTS ARE IN ONE CONTINUOUS B-REP. 
//     ( (K=4,N=5,T,RCOEF) IS THE TEMPORAL B-REP.) 
    taul[2] = tau[2];
    for(j=1; j<=ncd; ++j)
		rcoef[j*irc+3] = val[j*iv+1];
  
    for(i=1; i<=2; ++i){
		taul[i + 2] = tau[i + 2];
		for(j=1; j<=ncd; ++j)
			rcoef[i+3+j*irc] = val[i+1+j*iv];
		taul[3-i-1] = tau[1]-(tau[nm1]-tau[nm1-i]);
		for(j=1; j<=ncd; ++j)
		    rcoef[3-i+j*irc] = val[nv-i+j*iv];
    }
    *iflag=blg4sq_(3,3,taul,&rcoef[rcoef_offset],irc,5,ncd,irc,
		work, &t[1], &rcoef[rcoef_offset]);
    if(*iflag != 1)
		return;

// *** 3. OBTAIN TANGN(.). 
    ble_(4, 5, &t[1], &rcoef[rcoef_offset], irc, ncd, taul[2], 1, tangn);
// *** 4. INSERT TANGENT DATA AT THE BEGINING AND ENDING POINTS. 
    for(j=1; j<=ncd; ++j){
		rcoef[j*irc+1] = val[j*iv+1];
		rcoef[j*irc+2] = tangn[j-1];
		rcoef[nm1+j*irc] = tangn[j-1];
		rcoef[*n+j*irc] = val[nv+j*iv];
    }
    nvm1 = nv-1;
    for(i=2; i<=nvm1; ++i){
		for(j=1; j<=ncd; ++j)
			rcoef[i+1+j*irc] = val[i+j*iv];
    }
// *** 5. GET B-REP. 
    *iflag=blg4sq_(1,1,&tau[1],&rcoef[rcoef_offset],irc,*n,ncd,irc,
		work,&t[1],&rcoef[rcoef_offset]);
}
