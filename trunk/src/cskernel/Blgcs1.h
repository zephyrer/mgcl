/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLGCS1         BLGCS1 FOR BLGCS 
// BLGCS1 INSERTS CIRCLE DATA INTO GIVEN POINT SEQUENCE VALI(.,.), 
// AND GENERATE NEW KVAL(.,.) 
// *** INPUT *** 
// NCD.....SPACE DIMENSION OF VALI, MUST BE 2 OR 3 WHEN NK4>0. 
// NVI,KVALI(NVI),VALI(IVI,NCD)....INPUT POINT SEQUENCE WITH KNUCKLE 
//         INF.  NVI: NUM OF DATA            KVALI(.): KNUCKLE INF. 
//               VALI(.,NCD): POINTS OF NCD SPACE DIMENSION. 
// NK4,IDK(NK4),RCIR(NK4).....PROVIDE CIRCLE DATA. 
//       RCIR(I) IS RADIOUS OF THE OSCULATING CIRCLE AT THE IMAGINARY 
//       POINT VALI(L,.), I.E. KVALI(L)=1, 
//         WHERE L=IDK(I).       1<=I<=NK4 
//           IDK(.) MUST BE INCREASING, I.E. IDK(I)<IDK(I+1). 
// IVI,IV.....ROW DIMENSION OF VALI AND VAL, EACH. 
// *** OUTPUT *** 
// NV,KVAL(NV),VAL(IV,NCD)......OBTAINED NEW POINT SEQUENCE 
//         INCLUDING CIRCLE DATA. 
//         KVAL AND VAL MAY BE THE SAME AREA AS KVALI AND VALI, EACH. 
// IFLAG....=1: SUCCESSFUL RETURN     <>1:FAILURE 
// *** NOTE * 
//     WHEN KVAL(I)=-1, VAL(I,.) IS NOT POSITIONAL DATA, IS DERIVATIVE. 
//   BLGCS1 GENERATE -1  AS NEW KNUCKLE INF. THESE DATA ARE 
//   RESOLVED BY SUBROUTINE BLGCS3. SEE BLGCS3.
void blgcs1_(int ncd, int nvi, const int *kvali, 
	const double *vali, int nk4, const int *idk,const double *rcir, 
	int ivi, int iv, int *nv, int *kval, double *val, int *iflag
);
