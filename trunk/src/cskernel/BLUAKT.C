/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Blunk.h"
#include "cskernel/bk1fli.h"

// BLUAKT ADDS KNOT(S) INTO ORIGINAL B-REP, GUARANTREEING THE SAME LINE.
// *** INPUT  * 
//     K,N1,T1(N1+K),RCOEF1(IRC1,NCD),IRC1,NCD..... ORIGINAL B-REP. 
//     NAD,TAD(NAD),MLT(NAD).....NUM OF KNOT, PARAMETER VALUE, AND 
//              MULTIPLICITY AT THE PARAM TAD(.). 
//     IRC2....ROW DIMENSION OF RCOEF2 
// *** OUTPUT * 
//     N2,T2(N2+K),RCOEF2(IRC2,NCD)....NEW B-REP OBTAINED. 
// *** WORK   *   WORK(K,K) 
void bluakt_(int k, int n1, const double *t1, 
	const double *rcoef1, int irc1, int ncd, int nad, 
	const double *tad, const int *mlt, int irc2, double *work, 
	int *n2, double *t2, double *rcoef2
){
    int left, i, j, m, n1pk, leftp1;

    // Parameter adjustments 
    --t1;
    --mlt;
    --tad;
    --t2;

    // Function Body 
    n1pk = n1+k;
    for(i=1; i<=n1pk; ++i)
		t2[i] = t1[i];
    *n2 = n1;

    for(i=1; i<=nad; ++i){
		if(tad[i]>t2[k] && tad[i]<t2[*n2+1]){
			left=bk1fli_(*n2+1,&t2[1],tad[i]);
			m = mlt[i];
			j = left;
			while(t2[j]==tad[i]){
				--m; --j;
			}
			if(m>=1){
				leftp1 = left+1;
				for(j=*n2+k; j>=leftp1; --j)
					t2[j+m] = t2[j];
				for(j=1; j<=m; ++j)
					t2[left+j] = tad[i];
				*n2 += m;
			}
		}
    }
    blunk_(k,n1,&t1[1],rcoef1,irc1,ncd,*n2,&t2[1],irc2,work,rcoef2);
}
