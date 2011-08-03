/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//The original codes of this program comes from the FORTRAN code BANSLV of
//"A Practical Guide to Splines" by Carl de Boor.

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
void b1bslv_(const double *w, int nroww, int nrow, int nbandl, int nbandu, double *b){
    int jmax, i, j, nrowm1, middle;
    // Parameter adjustments 
    --b;
    w -= nroww + 1;

    // Function Body 
    middle = nbandu + 1;
	if(nrow>1){

    nrowm1 = nrow - 1;
    if(nbandl>0){
//                                 FORWARD PASS 
//            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
//            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) . 
    for (i = 1; i <= nrowm1; ++i) {
		jmax = nrow - i;
		if(jmax>nbandl)
			jmax = nbandl;
		for(j = 1; j <= jmax; ++j)
			b[i+j] -= b[i]*w[middle+j+ i*nroww];
    }
    }

//                                 BACKWARD PASS 
//            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
//            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN 
//            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW). 
    if (nbandu<=0){// A  IS LOWER TRIANGULAR . 
	    for(i=1; i<=nrow; ++i)
			b[i]/=w[i*nroww + 1];
	    return;
	}

	for(i=nrow; i>1; i--){
	    b[i] /= w[middle + i * nroww];
		jmax = i - 1;
		if(jmax>nbandu)
			jmax = nbandu;
		for (j = 1; j <= jmax; ++j)
			b[i-j]-=b[i]*w[middle-j+ i*nroww];
	}

	}
    b[1]/=w[middle+nroww];
}

