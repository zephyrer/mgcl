/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// bkcdtn_ WILL GENERATE DATA POINTS TAU(J) J=0,...,M-1 FROM 
// T(I) I=0,..N-1 SO AS TO BE PROPORTIONALLY DISTRIBUTED IN T. 
// *** INPUT *** 
// T(N) ..... OLD DATA POINT SEQUENCE (MAY BE KNOTVECTOR). 
// N ........ LENGTH OF T. 
// M ........ NUMBER OF NEW DATA POINTS TAU. 
// *** OUTPUT *** 
// TAU(M).... AREA TO STORE NEW DATA POINTS. 
void bkcdtn_(int n,const double *t, int m, double *tau){
    // Local variables 
    int nend;
    double tauc, dtau, span;
    int i, j, j1, j2;
    double tc, dt;

    // Parameter adjustments 
    --t; --tau;

    // Function Body 
    if (n <= 0 || m <= 0)
		return;
    tauc = t[1];
    if (t[n] <= t[1]) {
		for (i = 1; i <= m; ++i) tau[i] = tauc;
		return;
    }

// CREATE NEW DATA POINTS. 
    j1 = 1;
    j2 = n;
    if (m >= n) {
		while(j1 < n && t[j1] >= t[j1 + 1]) ++j1;
		if (j1 < n) {
		    while(j2 > 1 && t[j2 - 1] >= t[j2]) --j2;
		}
    }
    for (i = 1; i <= j1; ++i) tau[i] = tauc;
    nend = n - j2 + 1;
    for (i = 1; i <= nend; ++i) tau[m - i + 1] = t[n];

    span = t[j2] - t[j1];
    dt = span / (double) (j2-j1);
    if (m == 1) dtau = 0.f;
    else dtau = span / (double) (m+1-j1-nend);
	j1++;
    j = j1; j2 = m - nend;
    tc = tauc + dt;
    for (i = j1; i <= j2; ++i) {
		tauc += dtau;
		while(tc<=tauc) {tc+=dt; ++j;}
		tau[i] = t[j] - (t[j] - t[j - 1]) * ((tc - tauc) / dt);
    }
    return;
}
