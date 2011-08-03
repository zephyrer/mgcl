/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkktdp.h"
#include "cskernel/blgint.h"
#include "cskernel/blelin.h"
#include "cskernel/blg4sp2.h"

// BLUDKT GENERATES B-COEFFICIENTS OF DIFFERENT KNOT CONFIGURATION T2. 
// *** INPUT * 
//     K,N1,T1(N1+K),RCOE1(IRC1,NCD)IRC1,NCD...DESCRIBE THE OLD B-REP 
//          T1:KNOT VECTOR , RCOE1:B-COEF , N1:B-REP DIMENSION 
//     N2,T2(N2+K)....NEW BREP DIMENSION AND KNOT VECTOR 
//     IRC2....ROW DIMENSION OF THE VARIABLE RCOE2 
// *** OUTPUT * 
//     RCOE2(IRC2,NCD)..THE NEW B-COEF OBTAINED 
//     IFLAG... =1 :SUCCESFUL RETURN, <>1 : FAILURE BECAUSE OF ILLEGAL 
//                                          KNOT VECTOR T2. 
// *** WORK * 
//     WK1(N2),WK2(N2,9) 
// *** NOTE * 
//     THE NEW KNOT T2 MUST SATISFY THE FOLLOWLING CONDITION; 
//          T1(K)<=T2(K), AND T2(N2+1)<=T1(N1+1) 
void bludkt_(int k, int n1, const double *t1, 
	const double *rcoe1, int irc1, int ncd, int n2, 
	const double *t2, int irc2, double *wk1, double *wk2, 
	double *rcoe2, int *iflag)
{
    static int ibcbeg = 1;
    static int ibcend = 1;
    int rcoe2_offset, wk2_offset;
    int j;
    double tsave;

	// Parameter adjustments 
    wk2_offset = n2 + 1;
    wk2 -= wk2_offset;
    rcoe2_offset = irc2 + 1;
    rcoe2 -= rcoe2_offset;
// **************************** START OF BLUDKT  ************************
    bkktdp_(n2, k, t2, wk1);
    if(k==4){
		wk1[1] = wk1[0];
		wk1[n2-2] = wk1[n2-1];
		blelin_(k,n1,t1,rcoe1,irc1,ncd,ibcbeg,ibcend,
		n2,wk1,1,irc2,&rcoe2[rcoe2_offset]);
		for(j=1; j<=ncd; ++j){
		    tsave = rcoe2[j*irc2+1];
			rcoe2[j*irc2+1] = rcoe2[j*irc2+2];
		    rcoe2[j*irc2+2] = tsave;
		}
		*iflag = 2;
		for(j=1; j<=ncd; ++j){
			blg4sp2_(4,iflag,ibcbeg,ibcend,wk1,&rcoe2[j*irc2+1],
				irc2,n2,1,t2,1,&wk2[n2+1],&wk2[(n2<<1)+1],
				&wk2[n2*3+1],&rcoe2[j*irc2+1]);
		    if (*iflag!=1)
				break;
		}
		return;
    }else{
		blelin_(k,n1,t1,rcoe1,irc1,ncd,3,3, 
			n2,wk1,1,irc2,&rcoe2[rcoe2_offset]);
		*iflag=blgint_(wk1, &rcoe2[rcoe2_offset], t2, k, n2, ncd, irc2, irc2,
			 &wk2[wk2_offset], &rcoe2[rcoe2_offset]);
    }
}
