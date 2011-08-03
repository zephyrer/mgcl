/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkktdp.h"
#include "cskernel/blgint.h"
#include "cskernel/b1bslv.h"
#include "cskernel/bkcdtn.h"
#include "cskernel/bkmlt.h"
#include "cskernel/blelin.h"

// BLUDEC DECREASES B-REP DIMENSION ,APPPROXIMATING THE OLD B-REP. 
// *** INPUT * 
//     ISM   =1 This call of BUDEC is the first call as new data. 
//           =2 CALLED WITH SAME DATA as previous call. 
//              (Only RCOEF1 can be different from the previous call) 
//     K, N1,T1(N1+K),RCOEF1(IRC1,M),IRC1,M 
//            Original B-Spline. 
//     NDEC   : NUMBER OF KNOTS TO DECREASE 
//     IRC2   : ROW DIMENSION OF RCOEF2 
//     WORK1,WORK2 : WHEN ISM=2, PREVIOUS DATA MUST BE INPUT 
// *** OUTPUT * 
//     N2, T2(N2+K), RCOEF2(IRC2,M) 
//             Updated B-Spline. 
//     IFLAG  : = 1 GOOD 
//            : <> 1 ERROR 
// *** WORK   * 
//     WORK1(N2,2K-1) 
//     WORK2(N2) 
// *** NOTE *   . MAXIMUM MULTIPLICITY OF INNER KNOT IS K-1 
void bludec_(int ism, int k, int n1, 
	const double *t1, const double *rcoef1, int irc1, int m, 
	int ndec, int irc2, double *work1, double *work2, 
	int *n2, double *t2, double *rcoef2, int *iflag)
{
    int mltf, i;
	int nb, is;
    int nb2, km1, ndm, isf, ndz,nbm1;
    int rcoef2_offset, mltfm1;
    double rdec;

    // Parameter adjustments 
    --t1;
    rcoef2_offset = irc2 + 1;
    rcoef2 -= rcoef2_offset;
    --t2;

    // Function Body 
    km1 = k-1;
    if(ism!=2){
	
	//  --- CASE OF PROCESS FROM START 
	ndm = ndec;
	rdec = (double) (ndec) / (double) (n1 - k + 2);
	for(i=1; i<=k; ++i)
	    t2[i] = t1[i];
	*n2 = k;
	is = k;

	// ..LOOP UNTIL IS REACHES N1+1. 
	while(1){
		bkmlt_(n1,&t1[1],is,km1,&isf,&mltf);
		if(isf<0) 
		    isf=n1+1;
		nb = isf-is+1;
		ndz = (int)((double)nb*rdec);
		if(ndz>ndm)
		    ndz = ndm;
		nb2 = nb-ndz;
		ndm -= ndz;
		if(nb2<=1){
			++(*n2);
		    t2[*n2] = t1[isf];
		}else if(nb2!=nb){
			// DECREASE KNOT BY BKCDTN. 
			bkcdtn_(nb,&t1[is],nb2,&t2[*n2]);
		    *n2 = *n2+nb2-1;
		}else{
			// NO KNOTS DECREASED. 
			nbm1 = nb-1;
		    for (i = 1; i <= nbm1; ++i) {
			++(*n2);
			t2[*n2] = t1[is + i];
			}
		}

		if(isf>n1)
		    break;
		mltfm1 = mltf - 1;
		for(i=1; i<=mltfm1; ++i){
		    ++(*n2);
			t2[*n2] = t1[isf+i];
		}
		is = isf+mltfm1;
	}
	// ..END OF LOOP, ADD END KNOT 

	for (i = 1; i <= km1; ++i) {
	    ++(*n2);
	    t2[*n2] = t1[n1 + i];
	}
	*n2 -= k;
	//     NOW KNOT VECTOR GENERATED IN T2(.), 
	//     GENERATE DATA POINTS AND ASSOCIATED DATA IN (WORK2(.),RCOEF2()) 
	bkktdp_(*n2, k, &t2[1], work2);
	blelin_(k, n1, &t1[1], rcoef1, irc1, m, 3, 3, 
		*n2, work2, 1, irc2, &rcoef2[rcoef2_offset]);
	//      GET NEW B-REP. 
	*iflag=blgint_(work2, &rcoef2[rcoef2_offset], &t2[1], k, *n2, m, irc2, 
		irc2, work1, &rcoef2[rcoef2_offset]);

    }else{
	
	//   ---CASE OF USING LAST PROCESSED DATA. 
	blelin_(k, n1, &t1[1], rcoef1, irc1, m, 3, 3, 
		*n2, work2, 1, irc2, &rcoef2[rcoef2_offset]);
	for (i = 1; i <= m; ++i) {
	    b1bslv_(work1, k+km1, *n2, km1, km1, &rcoef2[i*irc2 + 1]);
	}
	*iflag = 1;
    
	}
}
