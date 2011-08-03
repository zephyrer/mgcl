/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/LinearEquation.h"
#include "mg/Matrix.h"
#include "mg/BPointSeq.h"

extern "C"{
#include "cskernel/b1bfac.h"
#include "cskernel/b1bslv.h"
#include "cskernel/b1hfac.h"
#include "cskernel/b1hslv.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// Solve linear equations.
//
//

//LU factorization using pivot.
//Factorize W to LU in the linear equation W*X=A.
//Using output W, the solution of W*X=A is obtained by solveLU.
//Here, L is n by n lower triangular matrix and U is n by n upper
//triangular matrix.  U's diagonal is (1).
void factorizeLU(
	MGMatrix& W,//left-hand side matrix of n*n is input, and
				//factorized LU matrix will be output where n=A.length().
				//left-bottom triangle including the diagonal will contain L,
				//right-upper triangle excluding the diagonal will contain U.
	int* id		//array of int id[n].
				//pivot id will be output in id[i] for i=0,...,n-1
){
	int i,j,k,n=W.sdim();
	int nm1=n-1, nm2=n-2;

	for(i=0; i<n; i++) id[i]=i;
	// LU factorization
	for(k=0; k<nm1; k++){
		int pivot_id=k;
		double pivot=W(k,k);
		for(j=k+1; j<n; j++){
			if(fabs(pivot)<fabs(W(j,k))){
				pivot_id=j; pivot=W(j,k);
			}
		}
		if(pivot_id!=k){
			for(j=0;j<n; j++){
				double w=W(k,j);
				W(k,j)=W(pivot_id,j);
				W(pivot_id,j)=w;
			}
			int iw=id[k];id[k]=id[pivot_id]; id[pivot_id]=iw;
		}
		for(j=k+1; j<n; j++) W(k,j)/=pivot;
		for(i=k+1;i<n; i++) for(j=k+1; j<n; j++) W(i,j)-=W(i,k)*W(k,j);
	}
}

//solve the linear equation W*X=A to get X inputting factorize W and A.
//W and id are obtained by factorizeLU.
void solveLU(
	const MGMatrix& W,
				//factorized LU matrix is input which is obtained by factorizeLU.
				//left-bottom triangle including the diagonal contain L,
				//right-upper triangle excluding the diagonal contain U.
	const int* id,	//array of int id[n] which contain pivot id in id[i] for i=0,...,n-1.
				//this is the output of factorizeLU.
	const MGBPointSeq& A,//right hand side vector. A.length() is W.sdim();
	MGBPointSeq& X	//solved X will be output. X.length() will be A.length(),
				//and X.sdim() will be A.sdim()
){
	int i,j,k,n=A.length();
	int m, sd=A.sdim();
	int nm1=n-1, nm2=n-2;

	X.resize(n,sd);
	for(m=0; m<sd; m++){//Loop for space dimension.
		for(i=0; i<n; i++) X(i,m)=A(id[i],m);
		//2. forward substitution
		for(k=0; k<n; k++){
			for(j=0;j<=k-1;j++) X(k,m)-=W(k,j)*X(j,m);
			X(k,m)/=W(k,k);
		}

		//3. backward substitution
		for(k=nm1; k>=0; k--){
			for(j=k+1; j<n; j++) X(k,m)-=W(k,j)*X(j,m);
		}
	}
}

//solve the linear equation W*X=A where W is a symetric tridiagonal matrix of order n.
void solveSymetricTridiagonal(
	MGBPointSeq& W,//the symetric tridiagonal matrix M is input.
		//three bands of the matrix at and above the diagonal are input.
		//W(i,0)=M(i,i),
		//W(i,1)=M(i,i+1) and M(i+1,i)
		//W(i,2)=M(i,i+2) and M(i+2,i) for i=0,...,n-1
		//Factorized matrix will be output.
	const MGBPointSeq& A,//right hand side vector. A.length() is W.length();
	MGBPointSeq& X	//solved X will be output. X.length() will be A.length(),
				//and X.sdim() will be A.sdim()
){
//	MGBPointSeq Wsave(W);
	int i,n=W.length();
	int j, sd=A.sdim();
	X.resize(n,sd);
	if(n<=1){
		for(j=0; j<sd; j++)//Loop for space dimension.
			X(0,j)=A(0,j)/W(0,0);
		return;
	}

//Solve W*X=A to get X.
	int nm1=n-1, nm2=n-2;

	//1. factorization
	for(i=0; i<nm1; i++){
		double ratio=W(i,1)/W(i,0);
		W(i+1,0)-=ratio*W(i,1);
		W(i+1,1)-=ratio*W(i,2);
		W(i,1)=ratio;
		ratio=W(i,2)/W(i,0);
		if(i<nm2) W(i+2,0)-=ratio*W(i,2);
		W(i,2)=ratio;
	}

	for(j=0; j<sd; j++){//Loop for space dimension.
		//2. forward substitution
		X(0,j)=A(0,j);
		X(1,j)=A(1,j)-W(0,1)*X(0,j);
		for(i=1; i<nm1; i++) X(i+1,j)=A(i+1,j)-W(i,1)*X(i,j)-W(i-1,2)*X(i-1,j);
		//3. backward substitution
		X(nm1,j)=X(nm1,j)/W(nm1,0);
		X(nm2,j)=X(nm2,j)/W(nm2,0)-X(nm1,j)*W(nm2,1);
		for(i=n-3; i>=0; i--) X(i,j)=X(i,j)/W(i,0)-X(i+1,j)*W(i,1)-X(i+2,j)*W(i,2);
	}
/*	cout<<endl<<"A="<<A<<"original_W*X=";
	for(i=0; i<n; i++){
		double wa=0.;
		int im2=i-2;if(im2>=0) wa+=X(im2,0)*Wsave(im2,2);
		int im1=i-1;if(im1>=0) wa+=X(im1,0)*Wsave(im1,1);
		wa+=X(i,0)*Wsave(i,0);
		int ip1=i+1;if(ip1<n) wa+=X(ip1,0)*Wsave(i,1);
		int ip2=i+2;if(ip2<n) wa+=X(ip2,0)*Wsave(i,2);
		cout<<i<<"="<<wa<<",";
	}
	cout<<endl;*/
}

/// factorizeBandLU executes the LU-factorization without pivoting of a banded matrix M.
/// M's order is n with (nlower + 1 + nupper) bands of diagonals stored in W.
/// GAUSS Elimination  W I T H O U T  pivoting is used. The routine is
/// intended for use with matrices M which do not require row 
/// interchanges durnig factorization, especially for the TOTALLY POSITIVE 
/// matrices which occur in spline calculations. 
/// The routine should not be used for arbitrary banded matrices. For arbitary ones, use
/// factorizeLU and solveLU.
/// On function's return, W contains the LU factorization of the matix M.
/// Function's return value is:
///              =0:SUCCESS, =1:FAILURE.
int factorizeBandLU(
	MGBPointSeq& W,///< contains interesting part of the matrix M.
///< W.length()=nlower+1+nupper=nband(band width of M), and 
///< W.sdim()=n(order of the matrix M).
///< The diagonals(or bands) of M(i+j,j) are stored in W as W(i+nupper,j)
///< for i=-nupper,..., nlower, and j=0,..., n-1.
///< Explicitly, M has nlower bands below the diagonal and nupper bands above the
///< diagonal. Thus the band width nband=nlower+1+nupper(=W.length()).
///< For exapmle, the interesting entries of M of order 9, whose nlower=1 and nupper=2
///< would appear in the 4 subscripts of W(i, j) as follows:
///<
///<     j=           0  1  2  3  4  5  6  7  8
///<     i=0:         x  x 02 13 24 35 46 57 68
///<     i=1:         x 01 12 23 34 45 56 67 78
///<     i=2:        00 11 22 33 44 55 66 77 88
///<     i=3:        10 21 32 43 54 65 76 87  x
///<
///< All other entries of W not identified in this way with an entry of M
///< are never referenced.
	int nlower	///<number of bands of the matrix M below the main diagonal.
){
	int n=W.sdim(), nband=W.length();
	assert(nband>nlower);
	assert(nlower>=0);
	return b1bfac_(W.data(), W.capacity(), n, nlower, nband-1-nlower);
}

/// solveBandLU is the companion routine of factorizeBandLU.
/// solveBandLU returns the solution of the linear system M*X = A
/// in place of M, given the LU-Factorization for M obtained by factorizeBandLU
/// in W.
///
/// ******  M E T H O D  ****** 
///   (With  M = L*U, as stored in  W,) the unit lower triangular system 
///  L(U*X) = B  is solved for Y = U*X, and  Y stored in B. Then the 
///  upper triangular syrtem  U*X = Y  is solved for X . The calculations 
///  are so arranged that the innermost loops stay within columns.
///  Here, B=A(.,i) and X=X(.,i) for i=0,...,A.sdim()-1.
void solveBandLU(
	const MGBPointSeq& W,///< factorized LU matrix is input
		///< which is obtained by factorizeBandLU.
		///< Left-bottom triangle including the diagonal contain L,
		///< right-upper triangle excluding the diagonal contain U.
		///< Refer to factorizeBandLU.
	int nlower,	///< number of bands of the matrix M below the main diagonal.
	const MGBPointSeq& A,///< right hand side vector. A.length() is W.sdim()
		///< =n(order of the matrix M).
	MGBPointSeq& X	///< solved X will be output. X.length() will be A.length(),
				///< and X.sdim() will be A.sdim().
){
	int n=W.sdim(), nband=W.length(), nsd=A.sdim(), nroww=W.capacity();
	assert(A.length()==n);
	int nupper=nband-1-nlower; assert(nupper>=0);
	const double* wd=W.data();
	X=A;
	for(int i=0; i<nsd; i++)
		b1bslv_(wd, nroww, n, nlower, nupper, &(X(0,i)));
}

/// Constructs Cholesky factorization.
/// Let C is a symmetric positive semidefinite and banded matrix, having nbands
/// diagonals at and below the main diagonal. Then C can be factorized as
///                   C = L*D*L-transpose
/// with L, unit lower triangular and D, the diagonal of C.
/// ******  M E T H O D  ****** 
/// Gauss elimination , adapted to the symmetry and bandedness of C, is used.
/// Near zero pivots are handled in a special way. The diagonal element
/// C(i,i) = W(0,i) is saved initially in wrok array diag(i), for all i.
/// At the i-th elimination step, the currrent pivot element, viz., W(0,i),
/// is compared with its original value, diag(i). If, as the result of prior
/// elimination steps, this element has been reduced by about a word length,
/// (i.e., if W(0,i)+dang(i) <= diag(i)), then the pivot is declared to be zero.
/// and the entire i-th row is declared to be linearly dependent on the prededing rows.
/// This has the effect of producing X(i)=0 when solving C*X=B for X, regardless of B.
/// Justification for this is as follows. In contemplated applications of this program,
/// the given equations are the normal equations for some least-squares approximation
/// problem, diag(i)=C(i,i) gives the norm-square of the i-th basis function, and,
/// at this point, W(0,i) contains the nrom-square of the error in the least-squares
/// approximation to the i-th basis function by linear combinations of the first i-1.
/// Having W(0,i)+diag(i) <= diag(i) signifies that the n-th function is
/// linearly dependent to machine accuracy on the first i-1 funcitons, therefore
/// can safely be left out from the basis of approximationg functions.
/// The solution of a linear system  C*X=B is effected by the succession of
/// the following two calls:
/// factorizeCholeLU to get factorization, and solveCholeLU  to solve for X.
void factorizeCholeLU(
	MGBPointSeq& W ///< contains nbands diagonals in W(.,j), with the main diagonal
///< in W(.,0). Precisely, W(i,j)  contains  C(i+j,j), i=0,...,nbands-1, j=0,...,n-1. 
///< W.length()=nbands, and W.sdim()=n is the order of the matrix C.
///< For example, the entries of a 7 diagonal symmetric matrix C of order 9
///< (nbands=4) would be stored in W(i,j) as
///<             j=          0  1  2  3  4  5  6  7  8
///<           i=0:         00 11 22 33 44 55 66 77 88   
///<           i=1:         10 21 32 43 54 65 76 87      
///<           i=2:         20 31 42 53 64 75 86         
///<           i=3:         30 41 52 63 74 85            
///< All other entries of W not identified in this way with an entry of C
///< are never referenced.
///< On return, W contains the Cholesky factorization C= L*D*L-transpose
///< with W(0,j) containing 1/D(j,j) for j=0,..., n-1.
///< And W(i,j) containing L(i+j,j) for i=1,..., nbands-1, and j=0, ..., n-i-1.
){
	int n=W.sdim(), nband=W.length(), ncapa=W.capacity();
	if(nband!=ncapa)
		W.reshape(nband);
	std::auto_ptr<double> work(new double[n]);
	b1hfac_(W.data(), nband, n, work.get());
}

/// Solves the linear system  C*X = A, provided W contains the Cholesky factorization.
/// The Cholesky factorization is obtained by factorizeCholeLU. See factorizeCholeLU. 
void solveCholeLU(
	const MGBPointSeq& W,///< contains the Cholesky factorization
		///< for C, can be obtained by factorizeCholeLU.
		///< Refer to factorizeCholeLU.
	const MGBPointSeq& A,///< right hand side vector. A.length() is W.sdim()=n(order of the matrix C).
	MGBPointSeq& X	///< solved X will be output. X.length() will be A.length(),
				///< and X.sdim() will be A.sdim().
){
	int n=W.sdim(), nband=W.length(), nsd=A.sdim();
	assert(A.length()==n);
	const double* wd=W.data();
	X=A;
	for(int i=0; i<nsd; i++)
		b1hslv_(wd, nband, n, &(X(0,i)));
}
