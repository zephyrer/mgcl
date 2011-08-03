/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bk1fli.h"

//bk2fli_ locates where input X is located in xt[.].
// *** INPUT * 
//   N,XT(N)....KNOT VECTOR OF LENGTH N. 
//   X        A POINT TO BE LOCATED IN XT 
// *** OUTPUT * 
//   LEFT...THE INDEX OF XT SUCH THAT 1<=LEFT<=N-1 and LEFT XT(LEFT) < X(LEFT+1), and
//          WHEN XT(1)<=X<=XT(N), XT(LEFT) <= X <= XT(LEFT+1). 
//          WHEN X<XT(1), LEFT is minimum integer s.t. X<XT(LEFT)<XT(LEFT+1).
//          WHEN XT(N)<X, LEFT is maximum integer s.t. XT(LEFT)<XT(LEFT+1)<X.
// *** NOTE *
//   x[0] must be less than x[n-1].
int bk2fli_(int n,const double *xt, double x){
    int lp1, left;
// ***************** START OF BK2FLI ********************************** 
	if(n<=1)
		return 1;

    // Parameter adjustments 
    --xt;

    left=bk1fli_(n, &xt[1], x);
    if(left==0){
		//FIND SMALLEST LEFT SUCH THAT XT(LEFT)<XT(LEFT+1)
		left = 1;
		while(left<n){
			lp1 = left + 1;
			if(xt[left]<xt[lp1])
				return left;
			left = lp1;
		}
    }else if(left>=n){
		//FIND LARGEST LEFT SUCH THAT XT(LEFT) < XT(LEFT+1)
		left = n - 1;
		lp1=n;
		while(1<left){
			if(xt[left]<xt[lp1])
				return left;
			lp1=left--;
		}
    }
	return left;
}
