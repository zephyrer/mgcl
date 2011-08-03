/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//b1bslv_ is the companion routine of b1bfac. It returns the solution of
//the linear system A*X = B in place of B, given the LU-Factorization for A in the array W.
//
// ******  I N P U T  ****** 
//  W, NROWW,NROW,NBANDL,NBANDU.....DESCRIBE THE LU-FACTORIZATION OF A 
//        BANDED MATRIX  A  OF RODER  NROW  AS CONSTRUCTED IN  B1BFAC . 
//        FOR DETAILS, SEE  B1BFAC . 
//  B.....RIGHT SIDE OF THE SYSTEM TO BE SOLVED . 
// ******  O U T P U T  ****** 
//  B.....CONTAINS THE SOLUTION  X , OF ORDER  NROW . 
// ******  M E T H O D  ****** 
//   (WITH  A = L*U, AS STORED IN  W,) THE UNIT LOWER TRIANGULAR SYSTEM 
//  L(U*X) = B  IS SOLVED FOR  Y = U*X, AND  Y  STORED IN  B . THEN THE 
//  UPPER TRIANGULAR SYSTEM  U*X = Y  IS SOLVED FOR  X  . THE CALCUL- 
//  ATIONS ARE SO ARRANGED THAT THE INNERMOST LOOPS STAY WITHIN COLUMNS.
void b1bslv_(const double *w, int nroww, int nrow, int nbandl, int nbandu, double *b);

