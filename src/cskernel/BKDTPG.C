/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//     bkdtpg_ WILL GENERATE DATA POINT SEQUENCE IN tau[I], 0 <=I<=n-1, 
//     S.T.          tau[0] = 0. 
//                   tau[I] = DIST(p(I,.)-p(I-1,.)) + tau[I-1] 
//                             FOR 1<=I<=n-1 . 
// *** INPUT * 
//     p(ip,ncd)     POINT-SEQUENCE IN A ncd-DIMENSION SPACE. 
//     n             NUM OF POINTS, p(I,.) 1<=I<=n  ARE EFFECTIVE. 
//     ncd           SPACE DIMENSION OF POINTS p( , ). 
//                   ( COLUMN DIMENSION OF THE VARIABLE p. ) 
//     ip            ROW DIMENSION OF THE VARIABLE p . 
// *** OUTPUT * 
//     tau(n)        DATA POINT SEQUENCE. 
void bkdtpg_(const double *p, int n, int ncd, int ip, double *tau){
    double sqrt(double);
    // Local variables 
    double d;
    int i, j;
    double d1;
    int im1;

    // Parameter adjustments 
    --tau;
    p -= ip+1;

    // Function Body 
    tau[1] = 0.f;
    for(i=2; i<=n; ++i){
		d = 0.f;
		im1 = i-1;
		for(j=1; j<=ncd; ++j){
		    d1 = p[i+j*ip] - p[im1+j*ip];
			d += d1*d1;
		}
		d = sqrt(d);
		tau[i] = tau[im1] + d;
    }
}
