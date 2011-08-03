/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLG4S1 IS A DEDICATED SUBROUTINE OF BLG4SQ. 
// BLG4S1 PRODUCES SHOENBERG'S VARIATION DIMINISHING KNOTS OF ORDER k.
// DATA POINTS TAU MAY BE MULTIPLE, I.E. 
// WHEN     TAU(I-1) < TAU(I)=---=TAU(I+M-1) < TAU(I+M)  , 
//          T(I+2+J) = TAU(I+J) FOR J=0,1,--,M-1 ( M<=3) . 
// *** INPUT * 
//     k..........order of the b-spline.
//     TAU(N).....DATA POINT SEQUENCE. 
//     N..........NUMBER OF DATA POINTS. 
// *** OUTPUT* 
//     T(N+4).....KNOT VECTOR OF SHOENBERG'S VARIATION DIMINISHING. 
//     IFLAG......INDICATES SUCCESS (=1) , OR FAILURE (=2). 
//          FAILURE IS CAUSED BY ILLEGAL NUMBER OF D.P. MULTIPLICITY; 
//           1) MORE THAN 3 MULTIPLICITY DETECTED. 
//           2) NO MUTIPLICITY ON STARTING (ENDING) AND NEXT (PREVIOUS) 
//              DATA POINT IS MULTIPLE. 
void blg4s1_(int k,const double *tau, int n, double *t, int *iflag){
    int i, j, jp1, jp2, jp3, jp4, nmk, km1;
	double tsum, fkm1;

    // Parameter adjustments 
    --t;
    --tau;

// **************************** START OF BLG4S1 ************************* 
    *iflag = 2;
    if(n<k)
		return;
    for(i=1; i<=k; ++i){
		t[i] = tau[1];
		t[n+i] = tau[n];
    }

// LOOP OVER J UNTIL J > N-K . 
    nmk = n-k;
    j = 1;
	if(k!=4){
		km1=k-1;
		fkm1=(double)km1;
		for(i=1; i<=nmk; i++){
			tsum=tau[i+1];
			for(j=2; j<=km1; j++)
				tsum+=tau[i+j];
			t[i+k]=tsum/fkm1;
		}
	    *iflag = 1;
		return;
	}

	while(j<=nmk){
	    jp1 = j+1;
		jp2 = j+2;
	    if(j>1){
			if(tau[jp1] >= tau[jp2])
				return;
		}
	    jp3 = j+3;
		jp4 = j+4;
		if(j>=nmk || tau[jp2]!=tau[jp3]){
			t[jp4] = (tau[jp1]+tau[jp2]+tau[jp3])/3.f;
			++j;
		    continue;
		}
		t[jp4] = tau[jp2];
		t[j+5] = tau[jp3];
		j += 2;
		if(tau[jp3]<tau[jp4])
			continue;
		t[j+4] = tau[j+2];
		++j;
	}
    *iflag = 1;
}
