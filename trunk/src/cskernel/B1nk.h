/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//   B1NK GETS K-COEF'S THAT ARE NECESSARY TO COMPUTE NEW J-TH B-COEF 
//  CORRESPONDING TO J-TH KNOT T(J). THE 'NEW' MEANS NEW SUBDIVIDED 
//  KNOT CONFIGURATION. 
// ** INPUT * 
//     K.........THE ORDER OF B-REP(NEW AND OLD). 
//     TAU(N+K)..THE OLD KNOT, WHERE N=B-REP DIMENSION OF THE OLD 
//     MU........THE SUBSCRIPT OF TAU S.T. TAU(MU) <= T(J) < TAU(MU+1). 
//                  (  MUST NOT  TAU(MU)=TAU(MU+1) ) 
//     T(M+K)....THE NEW KNOT, WHERE M=B-REP DIMENSION OF THE NEW 
//     J.........INDICATES AT WHICH COEF SHOULD BE EVALUATED AMONG 
//                    M B-COEF OF THE NEW B-REP. 
// ** OUTPUT * 
//     BATJ(K,K)=COEF'S THAT SHOULD BE MULTIPLIED TO RCOEF(.), 
//             I.E. BATJ(I,K) TO RCOEF(MU-K+I)  FOR  1<=I<=K. 
// ** WORK * 
//     BATJ(K,K)   WORK ARRAY OF LENGTH K*K( THE LAST K FOR OUTPUT AREA) 
// ** NOTE * 
// (1) B1NK EMPLOYS THE OSLO ALGORITHM OF COHEN,ET AL. 
// (2) J-TH B-COEF = SUM OF (BATJ(I,K)*RCOEF(MU-K+I)) 1<=I<=K, 
//     WHERE RCOEF IS  THE OLD B-COEF. 
void b1nk_(int k, const double *tau, int mu, const double *t, int j, double *batj);
