/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
int bk2fli_(int n,const double *xt, double x);
