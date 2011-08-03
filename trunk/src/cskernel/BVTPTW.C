/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Ble.h"
#include "cskernel/Bvtwst.h"

//  BVTPTW OBTAINS E(.,.), DEDU(.,.), AND DEDV(.,.) OF BVTWST FROM 
//  TANGENT PLANE (NTP,TTP,RTP) AND GETS TWIST DATA BY CALLING BVTWST. 
// ******  I N P U T  ****** 
//  LUD,UTAU(LUD),LVD,VTAU(LVD),SURF(ISR1,ISR2,3)........ 
//        POSITIONAL AND TANGENT DATA MISSING FOUR CORNER TWIST DATA. 
//        B-REP DIMENSION, DATA POINTS OF U AND V-DIRECTION EACH, 
//        AND DATA. 
//  IBC(4)....SPECIFIES BOUNDARY COND ALONG THE 4 BOUNDARIES, 
//        IBC(I)=1: TANGENT PLANE IN (NTP(I),TTP(ID(I,1),RTP(ID(I,2),.) 
//        IBC(I)<>1: NO T.P. 
//        ID I IS PERIMETER NUM. OF GIVEN PATCH. 
//  ID(4,2)...PROVIDE ID OF  TTP (ID(I,1)), AND 
//        RTP (ID(I,2)). THEY ARE VALID ONLY WHEN CORRESPONDING IBC(I) 
//        GIVES THE MEANING. 
//  NTP(4),TTP(.),RTP(IRTP,3).....PROVIDE TANGENT PLANE OF BOUNDARY, 
//        WHEN IBC(I)=1. 
//         (NTP(I),TTP(LT),RTP(LR,.)) IS THE TAN PLANE OF PERIMETER I, 
//             WHERE LT=ID(I,1), LR=ID(I,2). 
//  ISR1,ISR2,IRTP....ARE ARRAY LENGTH OF SURF (ISR1,ISR2), AND RTP 
// ***** W O R K ***** 
//  WK1(70),WK2(MM)  WORKARRAYS OF EACH LENGTH, WHERE MM=MAX(LUD,LVD) 
// *****  O U T P U T  ***** 
// SURF...CONTAINS TWIST DATA AT THEIR FOUR CORNER (SEE BVTWST) 
// IFLAG....=1 SUCCESSFUL RETURN,  <>1 SOME ERROR DETECTED. 
// ***** N O T E ***** 
// ORDER OF ALL B-REP ASSUMED TO BE 4. 
void bvtptw_(int lud,const double *utau, int lvd,const double *vtau,
	double *surf,const int *ibc,const int *id, 
	const int *ntp,const double *ttp,const double *rtp, int isr1,int isr2,
	int irtp, double *wk1, double *wk2, int *iflag)
{
    static int k = 4;
    double dedu[12],dedv[12],e[12];	//dedu,dedv,and e was [3][4].
    int i1, i2, nn;

    // Parameter adjustments 
    --ibc;
    id -= 5;
    --ntp;
    --ttp;
    rtp -= irtp+1;

//   <<<  OBTAIN E(.,.) , DEDU(.,.), AND DEDV(.,.) >>> 
    if(ibc[1]==1){
		i1 = id[5];
		i2 = id[9];
		nn = ntp[1];
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1],0,e);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1+nn],0,&e[3]);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1],1,dedu);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1+nn],1,&dedu[3]);
    }
    if(ibc[3]==1){
		i1 = id[7];
		i2 = id[11];
		nn = ntp[3];
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1],0,&e[9]);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1+nn],0,&e[6]);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1],1,&dedu[9]);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1+nn],1,&dedu[6]);
    }
    if(ibc[4]==1){
		i1 = id[8];
		i2 = id[12];
		nn = ntp[4];
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1],1,dedv);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1+nn],1,&dedv[9]);
    }
    if(ibc[2]==1){
		i1 = id[6];
		i2 = id[10];
		nn = ntp[2];
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1],1,&dedv[3]);
		ble_(k,nn,&ttp[i1],&rtp[i2+irtp],irtp,3,ttp[i1+nn],1,&dedv[6]);
    }
    bvtwst_(&ibc[1],e,dedu,dedv,lud,utau,lvd, vtau,surf,isr1,isr2,wk1,wk2,iflag);
}
