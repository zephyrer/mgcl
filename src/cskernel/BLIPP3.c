/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/blunk1.h"
#include "cskernel/blip2.h"
#include "cskernel/blip1.h"
#include "cskernel/bli1sp.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Store the parameter value tm in x[*nx].
//Function's return value is true if successfully stored. False if the area
//x[] overflowed.
int blipp31_store(double tm, double ts, int* nx, double* x, double errort, int mx){
    if(tm <= ts) return 1;
    if(*nx){
		if(tm - errort > x[*nx]) ++(*nx);
		if(*nx >= mx) return 0;
	}else{
		++(*nx);
	}
    x[*nx] = tm;
	return 1;
}

void blipp3_(int kin, int n,const double *t, 
	const double *rcoef, int *ki, double f, double error, 
	double errort, int mx, double ts, double *pstak, 
	double *tstak, double *work3, int *nx, double *x, 
	int *iend)
{
//         BLIPP3 WILL OBTAIN THE SOLUTION OF F=G(X) -- G:B-REP-- BY 
//         CONVERTING EACH SPAN INTO BEZIER CURVE AND SUBDIVIDING THE 
//         BEZIER CURVE. 
//         THE SOLUTION X (PARAM VALUE) ARE MORE THAN ONE,IN GENERAL. 
//         OSLO ALGORITHM IS EMPLOYED FOR SUBDIVISION. 
// *** INPUT * 
//     K,N,T(N+K),RCOEF(N).....PROVIDE B-REP. 
//     KI            INDEX OF RCOEF WHERE BLIPP3 TO START GETTING THE 
//                   SOLUTION. THE SOLUTION X IS >= T(KI+1) 
//     F             FUNCTION VALUE TO GET THE SOLUTION OF F=G(X) 
//     ERROR         ERROR ESTIMATE S.T. ABS(F-G(X)) <= ERROR 
//     ERRORT        ERROR ESTIMATE OF KNOT VECTOR. 
//     MX............MAX LENGTH OF TH EVARIABLE X(.). 
//     TS...........LOWER PARAMETER LIMIT, I.E. OBTAINED X(.)>TS 
//                  ( X(.) SUCH THAT X(.)<=TSIS DISCARDED. ) 
//     NX            NUMBER OF THE SOLUTIONS OBTAINED BEFORE THE CALL 
//                   OF BLIPP3. NUMBER OF NEW SOLUTIONS WILL BE ADDED 
//                   TO NX. 
// *** OUTPUT * 
//     NX,X(NX)..... THE SOLUTIONS ARE SET INTO X(I) AND NX UPDATED. 
//     KI............INDEX OF RCOEF IS RETURNED S.T. 
//                    RCOEF(KI-I) > F FOR I=0,--,K-2    , 
//                OR  RCOEF(KI-I) < F FOR I=0,--,K-2. 
//                   WHEN KI=N, COMPUTATION REACHED TO THE END. 
//     IEND..........INDICATES IF BLIPP3 PERFORMED FULL COMPUTATION 
//                  OR NOT. IEND=0 NOT, BECAUSE MX IS TOO SMALL. 
//                          IEND=1 PERFORMED FULL COMPUTATION 
// *** WORK * 
//     PSTAK(K,K+1)  WORK ARRAY OF LENGTH K*(K+1). 
//                   THE FIRST K*2 WORDS FOR TEMPORAL BEZIER CURVE. 
//                   THE LAST K*(K-1) WORDS FOR BEZIER CURVE STACK AREA. 
//     TSTAK(K,2,K+1) USED FOR KNOT VECTOR AREA OF THE PSTAK. 
//                   THE LAST K*(K-1) WORDS FOR BEZIER CURVE STACK AREA. 
//     WORK3(K*K) WORK ARRAY OF LENGTH K*K, FOR BLUNK1 AND BLI1SP. 
//                     USED FOR RCBNK1,BLI1SP. 

	int k=kin;
	int k2=k+k,k3=k*3,k4=k*4,k5=k*5,k6=k*6;
	int k3p1=k3+1;
    int kp1= k+1;
    int km1=k-1, km2=k-2;
//       KI: INDEX OF RCOEF  , IT: INDEX OT T. 
    int it= *ki + 1, itp1;
//ISTK: STACK POINTER, I.E. PSTAK(.,ISTK),TSTAK(.,.,ISTK) ARE THE TOP OF THE STACK.
    int istk=2, istkp1;
	int nspan,nints;
    int i;
	double t0,t1,tm;

	// Parameter adjustments 
    tstak -= k3p1;
    pstak -= kp1;
    --t;
    --x;

    *iend = 0;
    if (it < k) it = k;
    itp1 = it + 1;
//LOOP OVER KI UNTIL (K-2) SPANS AFTER RCOEF(KI) ARE GREATER OR LESS THAN F

	while(*ki < n){

//    GET NEXT INTERVAL OF THE ORIGINAL B-REP. 
//    CHECK THE INTERSECTION NUMBER OF THE NEXT (K-2) SPAN. 
    nspan = n - *ki;
    if(km2 && (nspan>km2)) nspan = km2;
    nints = blip1_(rcoef, *ki, nspan, f);
    if(nints == 0){
		*ki += nspan;
		*iend = 1;
		return;
	}

knot_index_loop://  LOOP OVER IT UNTIL  IT => KI+K. 
	while(it<=n && (t[it] >= t[itp1])){ it++;itp1++;}
    if(it>n) break;
	if(it>=(*ki+k)){ ++(*ki); continue;}
// CONVERT THE ORIGINAL B-REP TO BEZIER. 
	t0=t[it], t1=t[itp1];
    for (i = 1; i <= k; ++i) {
		tstak[i+k3] = t0;
		tstak[i+k4] = t1;
    }
    for(i=1; i<=k; ++i){
		pstak[i+k]=blunk1_(&t[1],rcoef,k,it,&tstak[k3p1],i,work3);
    }

// NOW PSTAK(.,1) & TSTAK(.,.,1) CONTAINS BEZIER CURVE, 
// GET THE INTERSECTION NUMBER OF THE BEZIER WITH F AND 
// BRANCH ACCORDING TO THE INTSCTN NO. 
stack_loop:	//Loop over stack(1 original knot span).
	t0=tstak[k4];
	t1=tstak[1+k4];
    tm = (t0+t1) * .5f;
    if(blip2_(&pstak[kp1],1,k,f) <= error){
	//  G(X)=F FOR TSTAK(K,1,1) <= X <= TSTAK(1,2,1). 
		if(!blipp31_store(tm,ts,nx,x,errort,mx)) return;
		goto pop_stack;
	}
    nints=blip1_(&pstak[kp1], 1, km1, f);
    if(nints == 0) goto pop_stack;
    if(nints == 1){
	// NOW THE SOLUTION IS THE ONLY ONE IN TSTKA(K,1,1)<=X<=TSTAK(1,2,1). 
	    tm=bli1sp_(k,k,&tstak[k3p1],&pstak[kp1],k,f,error,work3);
		if(!blipp31_store(tm,ts,nx,x,errort,mx)) return;
		goto pop_stack;
	};
	if((tm-t0)<errort){
		if(!blipp31_store(tm,ts,nx,x,errort,mx)) return;
		goto pop_stack;
	}
// NOW THE SOLUTION NUMBER MAY BE GRETER THAN 2, SUBDIVIDE THE BEZIER. 
    istkp1 = istk+1;
    for(i=1; i<=k; ++i){
		tstak[i+k5] = t0;
		tstak[i+k6] = tm;
		tstak[i+((istkp1<<1)+1)*k] = tm;
		tstak[i+((istkp1<<1)+2)*k] = t1;
    }
    for(i=1; i<=k; ++i){
		pstak[i+k2] =
			blunk1_(&tstak[k3p1],&pstak[kp1],k,k,&tstak[k5+1],i,work3);
		pstak[i+istkp1*k] =
			blunk1_(&tstak[k3p1],&pstak[kp1],k,k,&tstak[((istkp1<<1)+1)*k+1],i,work3);
    }
// SUBDIVIDED BEZIER CURVES ARE IN PSTAK(.,2) & PSTAK(.,ISTKP1) . 
    if(blip1_(&pstak[istkp1*k+1], 1, km1, f) >= 1) {istk = istkp1;}
//     TRANSFER THE FORMER SUBDIVIDED CURVE INTO PSTAK(.,1) . 
    for (i = 1; i <= k; ++i) {
		tstak[i+k3] = tstak[i+k5];
		tstak[i+k4] = tstak[i+k6];
		pstak[i+k] = pstak[i+k2];
    }
    goto stack_loop;
pop_stack:
//  GET THE BEZIER CURVE OF THE TOP OF THE STACK, AND GET THE SOLUTION. 
	if(istk<=2){
		it = itp1;
		itp1 = it + 1;
		goto knot_index_loop;
	}
//  POP THE STACKED BEZIER CURVE INTO PSTAK(.,1) . 
    for (i = 1; i <= k; ++i) {
		tstak[i+k3] = tstak[i+((istk<<1)+1)*k];
		tstak[i+k4] = tstak[i+((istk<<1)+2)*k];
		pstak[i+k] = pstak[i+istk*k];
    }
    --istk;
    goto stack_loop;
// CASE OF END OF THE SPAN. 

	}
    *ki = n;
    *iend = 1;
}

