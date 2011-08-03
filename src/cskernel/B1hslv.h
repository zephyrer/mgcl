/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  SOLVES THE LINEAR SYSTEM     C*X = B   OF ORDER  N R O W  FOR  X 
//  PROVIDED  W  CONTAINS THE CHOLESKY FACTORIZATION FOR THE BANDED (SYM- 
//  METRIC) POSITIVE DEFINITE MATRIX  C  AS CONSTRUCTED IN THE SUBROUTINE 
//    B 1 H F A C  (QUO VIDE). 
// ******  I N P U T  ****** 
//  NROW.....IS THE ORDER OF THE MATRIX  C . 
//  NBANDS.....INDICATES THE BANDWIDTH OF  C . 
//  W.....CONTAINS THE CHOLESKY FACTORIZATION FOR  C , AS OUTPUT FROM 
//        SUBROUTINE BCHFAC  (QUO VIDE). 
//  B.....THE VECTOR OF LENGTH  N R O W  CONTAINING THE RIGHT SIDE. 
// ******  O U T P U T  ****** 
//  B.....THE VECTOR OF LENGTH  N R O W  CONTAINING THE SOLUTION. 
// ******  M E T H O D  ****** 
//  WITH THE FACTORIZATION  C = L*D*L-TRANSPOSE  AVAILABLE, WHERE  L  IS 
//  UNIT LOWER TRIANGULAR AND  D  IS DIAGONAL, THE TRIANGULAR SYSTEM 
//  L*Y = B  IS SOLVED FOR  Y (FORWARD SUBSTITUTION), Y IS STORED IN  B, 
//  THE VECTOR  D**(-1)*Y IS COMPUTED AND STORED IN  B, THEN THE TRIANG- 
//  ULAR SYSTEM  L-TRANSPOSE*X = D**(-1)*Y IS SOLVED FOR  X (BACKSUBSTIT- 
//  UTION). 
void b1hslv_(const double *w, int nbands, int nrow, double *b);
