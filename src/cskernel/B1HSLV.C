/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

//The original codes of this program comes from the FORTRAN code BCHSLV of
//"A Practical Guide to Splines" by Carl de Boor.

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
void b1hslv_(const double *w, int nbands, int nrow, double *b)
{
    // Local variables 
    int jmax, j, n, nbndm1;

    // Parameter adjustments 
    --b;
    w -= nbands + 1;

    // Function Body 
    if(nrow<=1){
	    b[1] *= w[nbands + 1];
		return;
    }

//FORWARD SUBSTITUTION. SOLVE L*Y = B FOR Y, STORE IN B. 
    nbndm1 = nbands - 1;
    for(n=1; n<=nrow; ++n){
		jmax = nrow-n;
		if(jmax>nbndm1)
			jmax = nbndm1;
		for(j=1; j<=jmax; ++j)
			b[j+n] -= w[j+1+n*nbands]*b[n];
    }

//BACKSUBSTITUTION. SOLVE L-TRANSP.X = D**(-1)*Y  FOR X, STORE IN B. 
	for(n=nrow;n>0; --n){
	    b[n] *= w[n*nbands+1];
		jmax = nrow-n;
	    if(jmax>nbndm1)
			jmax = nbndm1;
		for(j=1; j<=jmax; ++j)
			b[n]-= w[j+1+n*nbands]*b[j+n];
	}
}
