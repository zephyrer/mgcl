/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"

// INTERNAL SUBROUTINE OF BLGCS, 
// BLGCS7 REMOVES DATA POINTS THAT ARE TOO CLOSE TO PREVIOUS DATA POINT 
// *** INPUT * 
//     TAU(N)...DATA POINT OF LENGTH N 
//     VAL(IV,NCD)..DATA POINT ORDINATE OF LENGTH N, 
//     KWORK(N).....INDICATES (TAU(I),VAL(I,.)) MAY BE REMOVED OR NOT: 
//          KWORK(I)=1 NOT ALLOWED TO REMOVE,  =0 ALLOWED. 
//     N........NUM OF DATA POINTS. 
//     NCD......SPACE DIMENSION OF THE DATA VAL 
//     IV.......ROW DIMENSION OF VAL 
// *** OUTPUT * 
//     TAU(N)....DATA POINTS UPDATED 
//     VAL(IV,NCD)...DATA POINT ORDINATES UPDATED 
//     N.....NEW NUM OF DATA POINTS 
void blgcs7_(double *tau, double *val, int *kwork, int *n, int ncd, int iv){
    // Local variables 
    int nnew, i, j;
    double dtold, ratio, dtnew;
    int nnewm1;
    if (*n <= 2)
		return;

    // Parameter adjustments 
    --kwork;
    --tau;
    val -= iv + 1;

    // Function Body 
    ratio = bkmax_();
    nnew = 2;
    i = 3;
	do{
	    if(kwork[nnew]==1){
			++nnew;
			kwork[nnew] = kwork[i];
			tau[nnew] = tau[i];
			for (j = 1; j <= ncd; ++j)
				val[nnew+j*iv] = val[i+j*iv];
		}else{
			nnewm1 = nnew-1;
			dtold = tau[nnew]-tau[nnewm1];
			dtnew = tau[i]-tau[nnew];
			if (dtold > dtnew * ratio) {
			//       CASE OF CURRENT SPAN IS TOO SHORT 
			    if (kwork[i] == 1) {
				//          DISCARD (TAU(NNEW),VAL(NNEW,.)). 
					kwork[nnew] = kwork[i];
					tau[nnew] = tau[i];
					for(j=1; j<=ncd; ++j)
						val[nnew+j*iv] = val[i+j*iv];
				}
			}else if(dtold*ratio<dtnew){
			//       CASE OF PREVIOUS SPAN IS TOO SHORT 
			    if(nnewm1>1 && kwork[nnewm1]!=1){
				//          DISCARD (TAU(NNEWM1),VAL(NNEWM1,.)) 
					kwork[nnewm1] = kwork[nnew];
					tau[nnewm1] = tau[nnew];
					for(j=1; j<=ncd; ++j)
						val[nnewm1+j*iv] = val[nnew+j*iv];
					nnew = nnewm1;
					--i;
			    }else{
				//          DISCARD (TAU(NNEW),VAL(NNEW,.)). 
					kwork[nnew] = kwork[i];
					tau[nnew] = tau[i];
					for(j=1; j<=ncd; ++j)
						val[nnew+j*iv] = val[i+j*iv];
				}
			}else{
				++nnew;
			    kwork[nnew] = kwork[i];
			    tau[nnew] = tau[i];
			    for(j=1; j<=ncd; ++j)
					val[nnew+j*iv] = val[i+j*iv];

			}
		}
	    ++i;
	}while(i <= *n);
    *n = nnew;
}
