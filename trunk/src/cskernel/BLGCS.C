/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/blg4sq.h"
#include "cskernel/blgcs1.h"
#include "cskernel/blgcs2.h"
#include "cskernel/blgcs3.h"

// BLGCS TO GEBERATE LINE B-REP ,GIVEN ITS POINT SEQUENCE WITH 
//  KNUCKLE INF. K. 
// *** INPUT * 
//       IBCI(2)......BOUNDARY COND. OF BOTH END POINTS, VALID ONLY 
//              WHEN KVAL(1) (OR KVAL(NV)) = 0. 
//                IBCI(.)=1: 1ST DERIV. PROVIDED IN TS(.) OR TE(.) 
//                       =2: 2ND DERIV. PROVIDED IN TS(.) OR TE(.) 
//                       =3: NO DERIV INF PROVIDED. 
//       TS(NCD),TE(NCD)....NECESSARY ONLY WHEN IBCI(.)=1 OR 2, AND 
//              KVAL(1) (KVAL(NV))=0, AND GIVES 1ST OR 2ND DERIVS. 
//       NCD....SPACE DIMENSION OF THE INPUT POINTS. 
//       NV....NUMBER OF INPUT POINTS, I.E. 
//              (  KVAL(I),VAL(I,.) )    FOR 1<=I<=NV 
//       KVAL(NV),VAL(IV,NCD)....INPUT POINTS WITH KNUCKLE INF. K 
//              ( KVAL(I),VAL(I,.) ) IS ONE PAIR OF THE INPUT POINTS 
//              KVAL(I) IS THE KNUCKLE INF. AND VAL(.,.) IS THE POINT 
//              (I.E. POSITIONAL DATA OF NCD SPACE DIMENSION) 
//       NK4,IDK(NK4),RCIR(NK4)....PROVIDE QSCULATING CIRCLE DATA 
//             NK4: NUM OF THE CIRCLE 
//             IDK(I):ID OF KVAL(.) WHERE KVAL(L)=3 AND RCIR(L)=RADIOUS 
//                    OF THE CIRCLE (L=KVAL(I)). 
//             RCIR(I):RADIOUS OF IDK(I). 
//       IV,IRC....ROW DIMENSIONS OF THE VARIABLE VAL AND RCOEF ,EACH. 
//                  IRC MUST BE .GE. OUTPUT N. 
// *** OUTPUT * 
//       N,T(N+4),RCOEF(IRC,NCD)....GENERATED LINE B-REP OF ORDER 4. 
//                RCOEF MAY BE THE SAME AREA AS VAL. 
//                N CAN BE APPROXIMATED BY; 
//                 N <= NV+2+NK4*7+(NUM OF KVAL(I)<>0) 
//       IFLAG...INDICATES IF SUCCESSFUL OR NOT, 
//                        = 1 SUCCESSFUL RETURN 
//                        <>1 FAILUE, T,RCOEF ARE INVALD. 
// ***WORK* 
//   WK1(N),WK2(MM)......WHERE MM=MAX(IRC*5+105,N*9) 
void blgcs_(const int *ibci, const double *ts, const double *te, int ncd,
	int nv,const int *kval,const double *val,int nk4,const int *idk,const double *rcir,
	int iv,	int irc, double *wk1, double *wk2, int *n, double *t, double *rcoef, int *iflag
){
    // System generated locals 
    int  wk2_offset;

    // Local variables 
    int ibc[2], nvo;

// **************************** START OF BLGCS *************************
// ===== 1. INSERT CIRCLE DATA ===== 
    // Parameter adjustments 
    wk2_offset = irc + 1;
    wk2 -= wk2_offset;

    // Function Body 
    blgcs1_(ncd,nv,kval,val,nk4,idk,rcir,iv, 
	    irc, &nvo, (int*)&wk2[irc+1], &wk2[(irc<<1)+1],iflag);
    if(*iflag == 1){
	// ===== 2. COMPUTE DATA POINT ===== 
		blgcs2_(ncd,nvo,(int*)&wk2[irc+1],&wk2[(irc<<1)+1],irc,&wk2[irc*5+1]);
		// ===== 3. INSERT DERIV INF OF STRAIGHT LINE ===== 
		blgcs3_(ibci,ts,te,ncd,nvo,&wk2[irc*5+1],(int*)&wk2[irc+1],&wk2[(irc<<1)+1],irc,irc,
			&wk2[irc*6+1],(int*)t,ibc,n,wk1,rcoef,iflag);
		if (*iflag == 1){
		// ===== 4. OBTAIN LINE B-REP ===== 
			*iflag=blg4sq_(ibc[0], ibc[1], wk1, rcoef, irc, *n, ncd, 
				irc, &wk2[wk2_offset], t, rcoef);
		}
    }
}
