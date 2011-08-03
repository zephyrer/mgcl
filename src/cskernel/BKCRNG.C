/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// bkcrng_ WILL CONVERT THE OLD KNOT SEQUENCE t[] INTO tnew[] SO THAT 
// tnew[k-1]=tstrt, tnew[n]=tend, AND MULTIPLICITY OF tnew[I] AND tnew[n+I] 
// HOLD. (0<=I<=k-1). 
//   WHEN k=0, t IS ASSUMED AS DATA POINT SEQUENCE. 
// *** INPUT * 
//     n,k,t[n+k]....INPUT KNOT VECTOR OF A B-REP, 
//                   n:B-REP DIMENSION, k:ORDER, t[.]:KNOT VECTOR 
//                   k MAY BE 0, AND INDICATES t[.] ARE DATA POINTS. 
//     tstrt,tend.......NEW PARAMETER RANGE OF NEW KNOT VECTOR,I.E. 
//                   VALUES SUCH THAT tnew[k-1]=tstrt AND tnew[n]=tend 
// *** OUTPUT * 
//     tnew[n+k].....NEW KNOT VECTOR OBTAINED 
//                   tnew MAY BE THE SAME AREA AS t. 
void bkcrng_(int k, int n, const double *t, double tstrt, double tend, double *tnew){
    // Local variables 
    double tatk,r;
    int i, i1, ie, is, np1, nmk;

    // Parameter adjustments 
    --tnew;
    --t;

    // Function Body 
    if(k>0){
		// CASE OF KNOT VECTOR CONVERT. 
		np1 = n+1;
	    tatk = t[k];
		r = (tend-tstrt) / (t[np1]-tatk);

	    nmk = n-k;
		for(i=1; i<=nmk; ++i)
			tnew[i+k] = tstrt+(t[i+k]-tatk)*r;
 
	//   COMPUTE MULTIPLICITY 
		for(is=1; is<k; ++is){
			if(tatk != t[k-is])
			    break;
	    }
	    for(ie=1; ie<k; ++ie) {
			if(t[np1] != t[np1+ie]){
			    break;
			}
		}

	    i1 =k-is;
	    for(i=1; i<=i1; ++i)
			tnew[i]= tstrt+(t[i]-tatk)*r;
	    for (i = 1; i <= is; ++i)
			tnew[i1 + i] = tstrt;

	    i1 = k-ie;
	    for(i=1; i<=i1; ++i)
			tnew[n+i+ie] = tend+(t[n+i+ie]-t[np1])*r;
		
	    for(i=1; i<=ie; ++i)
			tnew[n+i] = tend;

	}else{
	// CASE OF DATA POINT CONVERT 
	//   FIND MULTIPLICITY 
		for(is=2; is<=n; ++is){
			if(t[is] != t[1])
				break;
	    }
		ie = n-1;
		while(t[ie] == t[n])
			--ie;

	    r = (tend-tstrt)/(t[n]-t[1]);
	    if(is<=ie)
			for(i=is; i<=ie; ++i) 
				tnew[i] = tstrt+(t[i]-t[1])*r;

	    --is;
		++ie;
	    for(i=1; i<=is; ++i)
			tnew[i] = tstrt;
	    for(i=ie; i<=n; ++i)
			tnew[i] = tend;
	}
}
