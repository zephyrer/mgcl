/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BKKTDP WILL GENERATE DATA POINT SEQUENCE TAU(I), 0<=I<=N-1, FROM 
// KNOT VECTOR OF A B-REP. 
// *** INPUT * 
//     N,K,T(N+K)....INPUT KNOT VECTOR OF A B-REP, 
//                   N:B-REP DIMENSION, K:ORDER, T(.):KNOT VECTOR 
// *** OUTPUT * 
//     TAU(N).....DATA POINT SEQ. GENERATED. 
void bkktdp_(int n, int k, const double *t, double *tau){
    // Local variables 
    int i, j;
    double r;
    int id, km1, nm1, np1;
    double fkm1;

    // Parameter adjustments 
    --tau;
    --t;

    km1 = k-1;
    fkm1 = (double) km1;
    tau[1] = t[k];
    np1 = n+1;
    nm1 = n-1;
    for(i=2; i<=nm1; ++i){
		r = 0.f;
		for(j=1; j<=km1; ++j){
			id = i+j;
		    if(id<k)
				id = k;
		    if (id > np1)
				id = np1;
		    r += t[id];
		}
		tau[i] = r/fkm1;
    }
    tau[n] = t[np1];
}
