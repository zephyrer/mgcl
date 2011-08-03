/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bvdist.h"

// BLGCS2 is a dedicated subroutine of BLGCS, generates data point 
// sequence TAU(.) from data (VAL(j,NCD),j=1...N). VAL(j,.) may 
// include circle inf that is declared by KVAL(j). 
// ******Input****** 
//   NCD,N,KVAL(N),VAL(IV,NCD).....Input data of Space Dimension NCD, and 
//         of length N. KVAL(j) is a knuckle inf of VAL(j,.). 
// ******Output***** 
//   TAU(N).........Data point abssisa obtained of length N. 
void blgcs2_(int ncd, int n, const int *kval, 
	const double *val, int iv, double *tau
){
    // Local variables 
    double tnew;
    int i, j;
    double pntold[3], pntnew[3];

// **************************** START OF BLGCS2 *************************
    // Parameter adjustments 
    --tau;
    --kval;
    val -= iv+1;

    // Function Body 
    tnew = 0.f;
    tau[1] = tnew;
    for(j=1; j<=ncd; ++j)
		pntnew[j-1] = val[j*iv+1];

    for(i=2; i<=n; ++i){
		for(j=1; j<=ncd; ++j)
		    pntold[j-1] = pntnew[j-1];
		if(kval[i]==-1)
		    tau[i] = tnew;
		else{
		    for(j=1; j<=ncd; ++j)
				pntnew[j-1] = val[i+j*iv];
		    tnew += bvdist_(ncd, pntold, pntnew);
			tau[i] = tnew;
		}
    }
}
