/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Ble.h"
#include "cskernel/Bleval.h"

// BLELIN         BLELIN TO GENERATE INPUT DATA OF BLG4SQ FROM B-REP. 
// SUBROUTINE TO EVALUATE B-SPLINE(T(N+K),RCOEF(IRC,NCD)) 
// AT EACH DATA-POINT TAU(I),I=1,--,NT , AND CREATE VAL(NT,NCD) 
// *** INPUT * 
//  K,N,T(N+K),RCOEF(IRC,NCD),IRC,NCD....LINE B-REP TO EVALUATE 
//  IBCBEG.....INDICATES WHAT DATA BE GENERATED AT START 
//          =1 1ST DERIV    =2 2ND DERIV  =3 NO BC 
//  IBCEND.....INDICATES WHAT DATA BE GENERATED AT END 
//          =1 1ST DERIV    =2 2ND DERIV  =3 NO BC 
//  NT,TAU(NT)..DATA-POINTS 
//  IV1,IV2....1ST AND 2ND ARRAY LENGTH OF THE VARIABLE VAL. 
// *** OUTPUT * 
//  VAL(IV1,IV2,NCD).....VALUES EVALUATED 
// *** NOTE * 
//  EVALUATED VALUES ARE STORED IN ROW-WISE IN VAL, I.E. 
//         VAL(1,I,.) I=1....NT. 
void blelin_(int k, int n, const double *t, 
	const double *rcoef, int irc, int ncd, int ibcbeg, int ibcend,
	int nt,const double *tau, int iv1, int iv2, double *val
){
    // Local variables 
    int iend;
    double f[3];
    int i, j;
    int ntm1;

    // Parameter adjustments 
    --tau;
    val -= iv1*(iv2+1)+1;

    // Function Body 
    i = 1;
    if(ibcbeg!=3){
		ble_(k,n,t,rcoef,irc,ncd,tau[i],ibcbeg,f);
		for(j=1; j<=ncd; ++j)
			val[(i+j*iv2)*iv1+1] = f[j-1];
		++i;
    }
    ble_(k,n,t,rcoef,irc,ncd,tau[i],0,f);
    for(j=1; j<=ncd; ++j)
		val[(i+j*iv2)*iv1+1] = f[j-1];
    ++i;
    iend = nt-2;
    if(ibcend==3)
		iend = nt - 1;

	//   LOOP OVER I UNTIL I>IEND 
	while(i<=iend){
		//     ****     MULTIPLE POINT ? 
	    if (tau[i] != tau[i + 1]) {
		// NON MULTIPLE POINT. 
		//     ****     GET VALUES-OF-COORDINATES 
		    ble_(k, n, t, rcoef, irc, ncd, tau[i], 0, f);
			for(j=1; j<=ncd; ++j)
				val[(i+j*iv2)*iv1+1] = f[j-1];
			++i;
			continue;
	    }
		// CASE OF MULTIPLE DATA POINT (UP TO 3 MULTIPLICITY) 
		//    GET LEFT CONTINUOUS DERIVATIVE 
	    bleval_(k,n,t,rcoef,irc,ncd,1,tau[i],1,2,1,f);
		for(j=1; j<=ncd; ++j)
			val[(i+j*iv2)*iv1+1] = f[j-1];
		++i;
	    ble_(k,n,t,rcoef,irc,ncd,tau[i],0,f);
		for(j=1; j<=ncd; ++j)
			val[(i+j*iv2)*iv1+1] = f[j-1];
		++i;
	    if(tau[i]!=tau[i-1])
			continue;
		//    GET RIGHT CONTINUOUS DERIVATIVE 
	    ble_(k,n,t,rcoef,irc,ncd,tau[i],1,f);
	    for(j=1; j<=ncd; ++j)
			val[(i+j*iv2)*iv1+1] = f[j-1];
		++i;
	}
    if (ibcend != 3) {
	//      GET LEFT CONTINUOUS DERIVATIVE 
		ntm1 = nt-1;
		bleval_(k,n,t,rcoef,irc,ncd,1,tau[ntm1],ibcend,2,1,f);
		for(j=1; j<=ncd; ++j)
		    val[(ntm1+j*iv2)*iv1+1] = f[j-1];
    }
	bleval_(k,n,t,rcoef,irc,ncd,1,tau[nt],0,2,1,f);
    for(j=1; j<=ncd; ++j)
		val[(nt+j*iv2)*iv1+1] = f[j-1];
}
