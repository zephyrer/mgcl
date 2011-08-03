/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BSUDEC DECREASES NUMBER OF KNOTS OF SURFACE B-REP. 
// *** INPUT  * 
//      KU,LUD1,UKT1(LUD1+KU),KV,LVD1, 
//      VKT1(LVD1+KV),SURF1(IS11,IS12,3).....PROVIDE ORIGINAL SURF B-REP. 
//      KDEC,NDEC...INDICATE IN WHICH DIRECTION(KDEC) AND HOW MANY NUM 
//           (NDEC) SHOULD BE DECREMENTED. 
//           KDEC=1: U DIRECTION, =2: V. 
// *** OUTPUT * 
//      LUD2,UKT2(LUD2+KU),LVD2, 
//      VKT2(LVD2+KV),SURF2(IS21,IS22,3).....RETURN UPDATED SURF B-REP. 
//      IFLAG =1                     ..... NORMAL 
//            =2                     ..... ERROR (NEW. B. REP. CAL.) 
// *** WORK   * 
//      WORK1(MM,2*K-1),WORK2(MM),WORK3(LVD1,2), 
//             WHERE MM=MAX(LUD2,LVD2),K=MAX(KU,KV) 
void bsudec_(int ku, int lud1,const double *ukt1, 
	int kv, int lvd1,const double *vkt1,const double *surf1, 
	int kdec, int ndec, int is11, int is12, int is21, int is22,
	double *work1, double *work2, double *work3,
	int *lud2, double *ukt2, int *lvd2, double *vkt2, double *surf2, int *iflag
);
