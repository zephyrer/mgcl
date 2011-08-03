/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//b1bfac_ executes the LU-factorization without pivoting of the banded matrix A
//of order NROW with (NBANDL + 1 + NBANDU) bands of diagonals in W.
//  GAUSS Elimination  W I T H O U T  pivoting is used. The routine is
//  intended for use with matrices A which do not require row 
//  interchanges durnig factorzation, especially for the TOTALLY POSITIVE 
//  matrices which occur in spline calculations. 
//  The routine should not be used for arbitrary banded matrix.
//  Function's return value is:
//              =0:SUCCESS, =1:FAILURE.
//
// ******  I N P U T  ****** 
//  W.....ARRAY OF SIZE  (NROWW,NROW)  CONTAINING THE INTERESTING 
//      PART OF A BANDED MATRIX  A , WITH THE DIAGONALS OR BANDS OF  A 
//      STORED IN THE ROWS OF  W , WHILE COLUMNS OF  A  CORRESPOND TO 
//      COLUMNS OF  W . THIS IS THE STORAGE MODE USED IN  LINPACK  AND 
//      RESULTS IN EFFICIENT INNERMOST LOOPS. 
//         EXPLICITLY,  A  HAS  NBANDL  BANDS BELOW THE DIAGONAL 
//                          +     1     (MAIN) DIAGONAL 
//                          +   NBANDU  BANDS ABOVE THE DIAGONAL 
//      AND THUS, WITH    MIDDLE = NBANDU + 1, 
//        A(I+J,J)  IS IN  W(I+MIDDLE,J)  FOR I=-NBANDU,...,NBANDL 
//                                            J=1,...,NROW . 
//      FOR EXAMPLE, THE INTERESTING ENTRIES OF A (1,2)-BANDED MATRIX 
//      OF ORDER  9  WOULD APPEAR IN THE FIRST  1+1+2 = 4  ROWS OF  W 
//      AS FOLLOWS. 
//                        13 24 35 46 57 68 79   
//                     12 23 34 45 56 67 78 89   
//                  11 22 33 44 55 66 77 88 99   
//                  21 32 43 54 65 76 87 98      
//
//      ALL OTHER ENTRIES OF W  NOT IDENTIFIED IN THIS WAY WITH AN EN- 
//      TRY OF  A  ARE NEVER REFERENCED . 
//  NROWW.....ROW DIMENSION OF THE WORK ARRAY  W . 
//        MUST BE  .GE.  NBANDL + 1 + NBANDU  . 
//  NBANDL.....NUMBER OF BANDS OF  A  BELOW THE MAIN DIAGONAL 
//  NBANDU.....NUMBER OF BANDS OF  A  ABOVE THE MAIN DIAGONAL . 
//
// ******  O U T P U T  ****** 
//  W....CONTAINS THE LU-FACTORIZATION OF A INTO A UNIT LOWER TRIANGU- 
//       LAR MATRIX  L  AND AN UPPER TRIANGULAR MATRIX  U (BOTH BANDED) 
//       AND STORED IN CUSTOMARY FASHION OVER THE CORRESPONDING ENTRIES 
//       OF  A . THIS MAKES IT POSSIBLE TO SOLVE ANY PARTICULAR LINEAR 
//       SYSTEM  A*X = B  FOR  X  BY A 
//             CALL B1BSLV ( W, NROWW, NROW, NBANDL, NBANDU, B ) 
//       WITH THE SOLUTION X  CONTAINED IN  B  ON RETURN . 
//    IF function's return value is 1, then
//       ONE OF NROW-1, NBANDL,NBANDU FAILED TO BE NONNEGATIVE, OR ELSE 
//       ON OF THE POTENTIAL PIVOTS WAS FOUND TO BE ZERO INDICATING 
//       THAT  A  DOES NOT HAVE AN LU-FACTORIZATION. THIS IMPLIES THAT 
//       A  IS SINGULAR IN CASE IT IS TOTALLY POSITIVE . 
 int b1bfac_(double *w, int nroww, int nrow, int nbandl, int nbandu);
