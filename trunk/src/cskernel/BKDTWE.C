/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkdtpg.h"

//     BKDTWE WILL GENERATE DATA POINT SEQUENCE IN TAU(I), 1 <=I<=N, 
//     S.T.          TAU(1) = 0. 
//                   TAU(I) = DIST(P(I,.)-P(I-1,.)) + TAU(I-1) 
//                             FOR 2<=I<=N . 
//     IF IBCBEG,IBCEND <> 3 , MULTIPLE D.P. ADDED AT THE ASSOCIATED 
//     D.P. 
// *** INPUT * 
//     IBCBEG,IBCEND INDICATES BOUNDAY POINT INF ,I.E. 
//              =3:ALL DATA OF VAL ARE POSITIONAL DATA 
//             <>3:VAL(2,.) OR VAL(N-1,.) ARE DERIVATIVE INF, AND 
//                 MULTIPLE D.P. ADDED AT THE POINT 
//     VAL(IV,NCD)   POINT-SEQUENCE IN A NCD-DIMENSION SPACE. 
//     N             NUM OF POINTS, P(I,.) 1<=I<=N  ARE EFFECTIVE. 
//     NCD           SPACE DIMENSION OF POINTS P( , ). 
//                   ( COLUMN DIMENSION OF THE VARIABLE P. ) 
//     IV            ROW DIMENSION OF THE VARIABLE VAL 
// *** OUTPUT * 
//     TAU(N)        DATA POINT SEQUENCE. 
void bkdtwe_(int ibcbeg, int ibcend, const double *valin, int n, int ncd, int iv, double *tau)
{
    int j;
    double rsave[20];	// was [10][2]
    int i1, i2;
    int nm1;
	double* val=valin;//val will be updated after saved, then will be restored.

    // Parameter adjustments 
    val -= iv+1;;
    --tau;

    // Function Body 
    nm1 = n-1;
    i1 = 1;
    i2 = n;
    if(ibcbeg!=3){
	    for(j=1; j<=ncd; ++j){
			rsave[j-1] = val[j*iv+2];
			val[j*iv+2] = val[j*iv+1];
	    }
		i1 = 2;
	    for(j=1; j<=ncd; ++j){
			rsave[j+9] = val[nm1+j*iv];
			val[nm1+j*iv] = val[n+j*iv];
	    }
		i2 = nm1;
    }
    bkdtpg_(&val[i1+iv],i2-i1+1,ncd,iv,&tau[i1]);
    if (ibcbeg != 3) {
		tau[1] = tau[2];
	    for(j=1; j<=ncd;++j)
			val[j*iv+2] = rsave[j-1];
	    tau[n] = tau[nm1];
		for(j=1; j<=ncd; ++j)
			val[nm1+j*iv] = rsave[j+9];
    }
}
