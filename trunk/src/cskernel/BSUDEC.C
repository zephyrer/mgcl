/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bludec.h"

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
){
    int ludpku,lvdpkv;
    int i, j, l;
    int ism, irc1, irc2;

    // Parameter adjustments 
    work3 -= lvd1+1;;
    surf1 -= is11*(is12+1)+1;
    surf2 -= is21*(is22+1)+1;

    if(kdec==1){
	// *** Decrement dimension of U direction. *** 
		*lvd2 = lvd1;
		lvdpkv = lvd1+kv;
		for(i=0; i<lvdpkv; ++i)
			vkt2[i] = vkt1[i];
		ism = 1;
		irc1 = is11*is12;
		irc2 = is21*is22;
		for(i=1; i<=lvd1; ++i){
		    bludec_(ism,ku,lud1,ukt1,&surf1[(i+is12)*is11+1],irc1,3,ndec,
				irc2,work1,work2,lud2,ukt2,&surf2[(i+is22)*is21+1],iflag);
		    ism = 2;
		    if(*iflag==2)
				return;
		}
    }else{
	// *** Decrement dimension of V direction.*** 
		*lud2 = lud1;
		ludpku = lud1+ku;
		for(i=0; i<ludpku; ++i)
			ukt2[i] = ukt1[i];
		ism = 1;
		for(i=1; i<=*lud2; ++i){
		for(j=1; j<=3; ++j){
			for(l=1; l<=lvd1; ++l)
				work3[l+lvd1] = surf1[i+(l+j*is12)*is11];
			bludec_(ism,kv,lvd1,vkt1,&work3[lvd1+1],lvd1,1,ndec,
				lvd1,work1,work2,lvd2,vkt2,&work3[(lvd1<<1)+1],iflag);
			ism = 2;
			if (*iflag == 2)
				return;
			for(l=1; l<=*lvd2; ++l)
				surf2[i+(l+j*is22)*is21]=work3[l+(lvd1<<1)];
	    }
		}
    }
}
