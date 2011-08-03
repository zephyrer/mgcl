/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLINEAREQUATION_HH_
#define _MGLINEAREQUATION_HH_
/** @addtogroup ALGORITHM
 *  @{
 */

///
///

class MGMatrix;
class MGBPointSeq;

/** 
 *	@brief LU factorization to solve linear equations, general version.
 *  This is a general solver using pivotting, and so is inefficient
 *  compared with band matrix of diagonally dominant equations.
 *  If the equation is known as special equations as band matrix
 *  (as of diagonally dominant equations, or as Symetric Tridiagonal),
 *  use factorizeBandLU and solveBandLU, factorizeCholeLU and solveCholeLU,
 *  or solveSymetricTridiagonal.
 *
 *  factorizeLU factorizes W to LU in the linear equation W*X=A.
 *  Using output W, the solution of W*X=A is obtained by solveLU.
 *  Here, L is n by n lower triangular matrix and U is n by n upper
 *  triangular matrix.  U's diagonal is (1).
 *  LU factorization using pivot.
 */
void factorizeLU(
	MGMatrix& W,///<left-hand side matrix of n*n is input, and
				///<factorized LU matrix will be output where n=A.length(),
				///<left-bottom triangle including the diagonal will contain L,
				///<right-upper triangle excluding the diagonal will contain U.
	int* id		///<array of int id[n],
				///<pivot id will be output in id[i] for i=0,...,n-1
);

///Solve the linear equation W*X=A to get X, inputting factorize W and A.
///W and id are obtained by factorizeLU.
void solveLU(
	const MGMatrix& W,
				///<factorized LU matrix is input which is obtained by factorizeLU.
				///<left-bottom triangle including the diagonal contain L,
				///<right-upper triangle excluding the diagonal contain U.
	const int* id,	///<array of int id[n] which contain pivot id in id[i] for i=0,...,n-1,
				///<this is the output of factorizeLU.
	const MGBPointSeq& A,///<right hand side vector. A.length() is W.sdim().
	MGBPointSeq& X	///<solved X will be output. X.length() will be A.length(),
				///<and X.sdim() will be A.sdim().
);

///Solve the linear equation M*X=A where M is a symetric tridiagonal matrix of order n.
void solveSymetricTridiagonal(
	MGBPointSeq& W,///<the symetric tridiagonal matrix M is input.
		///<three bands of the matrix at and above the diagonal are input,
		///<W(i,0)=M(i,i),
		///<W(i,1)=M(i,i+1) and M(i+1,i),
		///<W(i,2)=M(i,i+2) and M(i+2,i) for i=0,...,n-1,
		///<Factorized matrix will be output.
	const MGBPointSeq& A,///<right hand side vector. A.length() is W.length().
	MGBPointSeq& X	///<solved X will be output. X.length() will be A.length(),
				///<and X.sdim() will be A.sdim().
);

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
///
/// NOTE: This program is translated to C++ from BANFAC of the book
/// "A Practical Guide to Splines" by Carl de Boor, Springer-Verlag.
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
///<     j=           0  1  2  3  4  5  6  7  8:
///<     i=0:         x  x 02 13 24 35 46 57 68,
///<     i=1:         x 01 12 23 34 45 56 67 78,
///<     i=2:        00 11 22 33 44 55 66 77 88,
///<     i=3:        10 21 32 43 54 65 76 87  x.
///<
///< All other entries of W not identified in this way with an entry of M
///< are never referenced.
	int nlower	///<number of bands of the matrix M below the main diagonal.
);

/// solveBandLU returns the solution of the linear system M*X = A.
/// solveBandLU is the companion routine of factorizeBandLU, and
/// factorizeBandLU executes the LU-Factorization for M in W.
///
/// NOTE: This program is translated to C++ from BANSLV of the book
/// "A Practical Guide to Splines" by Carl de Boor, Springer-Verlag.
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
);

/// Constructs Cholesky factorization.
/// Let C is a symmetric positive semidefinite and banded matrix, having nbands
/// diagonals at and below the main diagonal. Then C can be factorized as
///                   C = L*D*L-transpose
/// with L, unit lower triangular and D, the diagonal of C.
///
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
///
/// NOTE: This program is translated to C++ from BCHFAC of the book
/// "A Practical Guide to Splines" by Carl de Boor, Springer-Verlag.
void factorizeCholeLU(
	MGBPointSeq& W ///< contains nbands diagonals in W(.,j), with the main diagonal
///< in W(.,0). Precisely, W(i,j) contains  C(i+j,j), i=0,...,nbands-1, j=0,...,n-1. 
///< W.length()=nbands, and W.sdim()=n is the order of the matrix C.
///< For example, the entries of a 7 diagonal symmetric matrix C of order 9
///< (nbands=4) would be stored in W(i,j) as follows:
///<
///<             j=          0  1  2  3  4  5  6  7  8:
///<           i=0:         00 11 22 33 44 55 66 77 88,
///<           i=1:         10 21 32 43 54 65 76 87  x,
///<           i=2:         20 31 42 53 64 75 86  x  x,       
///<           i=3:         30 41 52 63 74 85  x  x  x.          
///<
///< All other entries of W not identified in this way with an entry of C
///< are never referenced.
///< On return, W contains the Cholesky factorization C= L*D*L-transpose
///< with W(0,j) containing 1/D(j,j) for j=0,..., n-1.
///< And W(i,j) containing L(i+j,j) for i=1,..., nbands-1, and j=0, ..., n-i-1.
);

/// Solves the linear system  C*X = A , provided W contains the Cholesky
/// factorization.
/// The Cholesky factorization is obtained by factorizeCholeLU. See factorizeCholeLU. 
///
/// NOTE: This program is translated to C++ from BCHSLV of the book
/// "A Practical Guide to Splines" by Carl de Boor, Springer-Verlag.
void solveCholeLU(
	const MGBPointSeq& W,///< contains the Cholesky factorization
		///< for C, can be obtained by factorizeCholeLU.
		///< Refer to factorizeCholeLU.
	const MGBPointSeq& A,///< right hand side vector. A.length() is W.sdim()=n(order of the matrix C).
	MGBPointSeq& X	///< solved X will be output. X.length() will be A.length(),
				///< and X.sdim() will be A.sdim().
);

/** @} */ // end of ALGORITHM group
#endif
