/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bk1fli.h"

// BKDMIX WILL MIX DATA SEQ. T1 AND T2, STORED IN T. 
// *** INPUT * 
//  N1          NUM OF DATA SEQ T1. 
//  T1(N1)      DATA SEQ 1. 
//  N2          NUM OF DATA SEQ T2. 
//  T2(N2)      DATA SEQ 2. 
// *** OUTPUT * 
//  N           NUM OF NEW DATA SEQ T. 
//  T(N)        NEW DATA SEQUENCE. 
// *** NOTE * 
//    DATA POINT MULTIPLICITY IS ALLOWED ONLY IN T1. 
void bkdmix_(int n1, const double *t1, int n2, const double *t2, unsigned *n, double *t){
    // Local variables 
    int j2mj1, i, i1, j1, j2, k2, k1, i2, istrt;
    int ii, im, j2m, j2m1, j1p1, j1p2;

    // Parameter adjustments 
    --t1;
    --t2;
    --t;

    // Function Body 
    j1=bk1fli_(n2, &t2[1], t1[1]);
    if (j1 < 1)
		j1 = 1;
    *n = 1;
    t[*n] = t1[1];
    i = 2;

// LOOP OVER I UNTIL I GT N1 
    while(i<=n1){
		j2=bk1fli_(n2, &t2[1], t1[i]);
		if(j2<1)
		    j2=1;
		if(t1[i] == t2[j2]){
			k2 = 1;
		    j2m = j2;
		    while(j2m >= 2 && t2[j2m] == t2[j2m - 1]) {
			++k2;
			--j2m;
		    }
			k1 = 1;
		    im = i;
		    while(im >= 2 && t1[im] == t1[im - 1]) {
			++k1;
			--im;
			}
		    if(k2-k1>0)
				j2 -= k2 - k1;
		}
		j2mj1 = j2-j1;
		if(j2mj1>1){
//         CASE OF T1(I) JUMPED MORE THAN 1 DATA POINTS. 
		    j1p2 = j1+2;
		    j1p1 = j1+1;
		    j2m1 = j2-1;
			if (j2mj1 <= 2) {
				//CASE OF J2-J1=2. 
				++(*n);
				t[*n] = (t[*n-1]+t2[j1p1]+t2[j1p2]+t1[i])/4.f;
		    }else{
				//CASE OF J2-J1 >= 3. 
				if(t2[j1p1]-t[*n] > t1[i]-t2[j2]){
				    istrt = 1;
				    ++(*n);
					t[*n] = (t[*n-1]+t2[j1p1]+t2[j1p2])/3.f;
					i1 = j1p2;
				    i2 = j2 - 2;
				}else{
				    istrt = 0;
				    ++(*n);
				    t[*n] = (t[*n-1]+t2[j1p1]+t2[j1p2]+t2[j1+3])/4.f;
					i1 = j1 + 3;
					i2 = j2m1;
				}
			//INSERT T2 DATA OF BETWEEN 
				if (i1 <= i2) {
				    for (ii = i1; ii <= i2; ++ii) {
					++(*n);
					t[*n] = t2[ii];
					}
				}
				if(istrt == 1){
				    ++(*n);
				    t[*n] = (t2[j2-2]+t2[j2m1]+t2[j2]+t1[i])/4.f;
				}else{
				    ++(*n);
				    t[*n] = (t2[j2m1]+t2[j2]+t1[i])/3.f;
				}
			}
		}

		++(*n);
		t[*n] = t1[i];
		++i;
		j1 = j2;
    }
}
