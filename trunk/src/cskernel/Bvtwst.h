/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  BVTWST COMPUTES TWIST DATA FROM SURFACE DATA SURF(.,.,.), AND 
//  STORES BACK IN SURF. 
// ******  I N P U T  ****** 
//  IBC(4)....SPECIFIES BOUNDARY COND ALONG THE 4 BOUNDARIES, 
//        IBC(I)=1: TANGENT VECTOR ALONG THE PERIMETER I. 
//        IBC(I)<>1: NO B.C. FOR THE PERIMETER I. 
//        ID I IS PERIMETER NUM. OF GIVEN PATCH. 
//        WHEN IBC(I) <> 1, FOLLOWING E(.,I), DEDU(.,I), DEDV(.,I) 
//        ARE INVALID. 
//  E(3,4).....ARE NORMAL VECTORS OF THE SURFACE AT THE FOUR CORNER OF 
//        THE SURFACE. 
//  DEDU(3,4), DEDV(3,4)....ARE DERIVATIVES OF E(,.) ALONG U(DEDU) 
//         AND ALONG V (DEDV). 
//         I=1: U-MIN AND V-MIN, I=2: U-MAX AND V-MIN, 
//          =3: U-MAX AND V-MAX,  =3: U-MIN AND V-MAX 
//  LUD,UTAU(LUD),LVD,VTAU(LVD),SURF(ISR1,ISR2,3)....... 
//        TANGENT DATA INPUT AS FOLLOWS, ACCORDING TO IBC(I); 
//     SURF(I,2,.)     1<=I<=LUD ALONG V=MIN PARAM LINE(TAN VEC OF V DIR). 
//     SURF(LUD-1,J,.) 1<=J<=LVD ALONG U=MAX PARAM LINE(TAN VEC OF U DIR). 
//     SURF(I,LVD-1,.) 1<=I<=LUD ALONG V=MAX PARAM LINE(TAN VEC OF V DIR). 
//     SURF(2,J,.)     1<=J<=LVD ALONG U=MIN PARAM LINE(TAN VEC OF U DIR). 
//  ISR1,ISR2....ARE ARRAY LENGTH OF SURF (ISR1,ISR2) 
// ***** W O R K ***** 
//  WK1(70),WK2(MM)  WORKARRAYS OF EACH LENGTH, WHERE MM=MAX(LUD,LVD) 
// *****  O U T P U T  ***** 
// SURF...CONTAINS TWIST DATA AT THEIR FOUR CORNER; 
//     SURF(2,2,.)          AT U=MIN,V=MIN. 
//     SURF(LUD-1,2,.)      AT U=MAX,V=MIN. 
//     SURF(LUD-1,LVD-1,.)  AT U=MAX,V=MAX. 
//     SURF(2,LVD-1,.)      AT U=MIN,V=MAX. 
// IFLAG....=1 SUCCESSFUL RETURN, 
//          =85 SET AS TWIST DATA = ALL 0.0 
//          <>1 SOME ERROR DETECTED. 
// ***** N O T E ***** 
// ORDER OF ALL B-REP ASSUMED TO BE 4. 
void bvtwst_(const int *ibc,const  double *e,const double *dedu, 
	const double *dedv, int lud,const double *utau, int lvd,const double *vtau,
	double *surf, int isr1, int isr2, 
	double *wk1, double *wk2, int *iflag);
