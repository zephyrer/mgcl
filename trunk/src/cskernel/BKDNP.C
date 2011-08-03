/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BKDNP deletes too near data point(s). 
//
// ***** Input ****** 
//  N,TAU(N), VAL(IV,NCD) 
//        Data points of length N, 
//        (TAU(i),VAL(i,j)) is the data point sequence, must be 
//        non-decreasing. TAU(i) is 
//        the abssisa, and VAL(i,j) is the ordinates, 1 <= i <= N. 
//        NCD is the space dimension, and 1 <= j <= NCD. 
//        IV is the row dimension of the variable VAL. 
// IMLT   Indicates whether multiple data points be deleted or not. 
//        Multple data points mean the points, 
//           TAU(i)=TAU(i+1). 
//        IMLT =1  : delete multiple points. 
//             <>1 : not delete multiple data points. 
// RATIO  Ratio to regard as a too near point. ratio>1.
//        Let r = (TAU(i)-TAU(i-1))/(TAU(i-1)-TAU(i-2)), then 
//        if r > RATIO, TAU(i-1) is removed. 
//        if r < 1/RATIO, TAU(i-1) is removed.
// Let dtold=tau[i-1]-tau[i-2], and dtnew=tau[i]-tau[i-1], then
// if dnew>ratio*dtold, tau[i-2] will be removed.
// if dtold>dtnew*ratio, tau[i] will be removed.
//
// ***** Output ***** 
// N,TAU(N), VAL(IV,NCD)   will be updated. 
//
// ***** NOTE ***** 
// Start and End Points are never removed, instead the next point is 
// removed. 
void bkdnp_(int *n, double *tau, double *val, int iv, int ncd, int imlt, double ratio){
    int nnew, nold, i, j;
    double dtold, dtnew;
    int nnewm1;

	nold=*n;
    if(nold<=2)
		return;

    // Parameter adjustments 
    --tau;
    val -= iv+1;

    nnew = 2;
    i = 3;
    if(imlt==1){

// Case of multiple points removal, i.e. discard same points. 
	while(i<=nold){
		nnewm1 = nnew-1;
		dtold = tau[nnew]-tau[nnewm1];
		dtnew = tau[i]-tau[nnew];
		if(dtnew==0. || dtold>dtnew*ratio){
		//CASE OF CURRENT SPAN IS TOO SHORT. 
			if(i<nold){
				i++;
				continue;//tau[i] is discarded if i !=nold.
			}else{
				//discard tau[nold-1] instead of tau[i] since i==nold.
				tau[nnew] = tau[i];
				for(j=1; j<=ncd; ++j)
					val[nnew+j*iv] = val[i+j*iv];
				break;
			}
		}
		
		if(dtnew>dtold*ratio){
		//CASE OF PREVIOUS SPAN IS TOO SHORT, or too large current span is encountered.
			if(nnewm1<=1){
			// DISCARD (TAU(NNEW),VAL(NNEW,.)) since nnewm1==1.
				tau[nnew] = tau[i];
				for(j=1; j<=ncd; ++j)
					val[nnew+j*iv] = val[i+j*iv];
				i++;
				continue;;
			}else{
				// DISCARD (TAU(NNEWM1),VAL(NNEWM1,.)). 
				tau[nnewm1] = tau[nnew];
				for (j = 1; j <= ncd; ++j) {
					val[nnewm1+j*iv] = val[nnew+j*iv];
				}
				nnew = nnewm1;
				continue;
			}
		}

		//not short span. Check next tau.
		++nnew;
		tau[nnew] = tau[i];
		for(j=1; j<=ncd; ++j)
			val[nnew+j*iv] = val[i+j*iv];
		i++;
	}

    }else{

// Case of multiple points non-removal. 
	while(i<=nold){
		nnewm1 = nnew - 1;
		if(tau[nnew] != tau[nnewm1] && tau[i] != tau[nnew]){

		dtold = tau[nnew]-tau[nnewm1];
		dtnew = tau[i] - tau[nnew];

		if(dtold>dtnew*ratio){
		// CASE OF CURRENT SPAN IS TOO SHORT. 
			if(i==nold || tau[i] == tau[i+1]){
			//DISCARD (TAU(NNEW),VAL(NNEW,.)). 
				tau[nnew] = tau[i];
				for(j=1; j<=ncd; ++j)
					val[nnew+j*iv] = val[i+j*iv];
			}
			i++;
			continue;
		}

		if(dtnew>dtold*ratio){
		//CASE OF PREVIOUS SPAN IS TOO SHORT. 
			if(nnewm1<=1 || tau[nnewm1]==tau[nnewm1-1]){
			//DISCARD (TAU(NNEW),VAL(NNEW,.)). 
				tau[nnew] = tau[i];
				for(j=1; j<=ncd; ++j)
					val[nnew+j*iv] = val[i+j*iv];
				i++;
			}else{
			//DISCARD (TAU(NNEWM1),VAL(NNEWM1,.)). 
				tau[nnewm1] = tau[nnew];
				for(j=1; j<=ncd; ++j)
					val[nnewm1+j*iv] = val[nnew+j*iv];
				nnew = nnewm1;
			}
			continue;
		}

		}

		++nnew;
		tau[nnew] = tau[i];
		for(j=1; j<=ncd; ++j)
			val[nnew+j*iv] = val[i+j*iv];
		i++;
	}

	}

	*n = nnew;
	return;
}
