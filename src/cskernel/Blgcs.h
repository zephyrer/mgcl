/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLGCS TO GEBERATE LINE B-REP ,GIVEN ITS POINT SEQUENCE WITH 
//  KNUCKLE INF. K. 
// *** INPUT * 
//       IBCI(2)......BOUNDARY COND. OF BOTH END POINTS, VALID ONLY 
//              WHEN KVAL(1) (OR KVAL(NV)) = 0. 
//                IBCI(.)=1: 1ST DERIV. PROVIDED IN TS(.) OR TE(.) 
//                       =2: 2ND DERIV. PROVIDED IN TS(.) OR TE(.) 
//                       =3: NO DERIV INF PROVIDED. 
//       TS(NCD),TE(NCD)....NECESSARY ONLY WHEN IBCI(.)=1 OR 2, AND 
//              KVAL(1) (KVAL(NV))=0, AND GIVES 1ST OR 2ND DERIVS. 
//       NCD....SPACE DIMENSION OF THE INPUT POINTS. 
//       NV....NUMBER OF INPUT POINTS, I.E. 
//              (  KVAL(I),VAL(I,.) )    FOR 1<=I<=NV 
//       KVAL(NV),VAL(IV,NCD)....INPUT POINTS WITH KNUCKLE INF. K 
//              ( KVAL(I),VAL(I,.) ) IS ONE PAIR OF THE INPUT POINTS 
//              KVAL(I) IS THE KNUCKLE INF. AND VAL(.,.) IS THE POINT 
//              (I.E. POSITIONAL DATA OF NCD SPACE DIMENSION) 
//       NK4,IDK(NK4),RCIR(NK4)....PROVIDE QSCULATING CIRCLE DATA 
//             NK4: NUM OF THE CIRCLE 
//             IDK(I):ID OF KVAL(.) WHERE KVAL(L)=3 AND RCIR(L)=RADIOUS 
//                    OF THE CIRCLE (L=KVAL(I)). 
//             RCIR(I):RADIOUS OF IDK(I). 
//       IV,IRC....ROW DIMENSIONS OF THE VARIABLE VAL AND RCOEF ,EACH. 
//                  IRC MUST BE .GE. OUTPUT N. 
// *** OUTPUT * 
//       N,T(N+4),RCOEF(IRC,NCD)....GENERATED LINE B-REP OF ORDER 4. 
//                RCOEF MAY BE THE SAME AREA AS VAL. 
//                N CAN BE APPROXIMATED BY; 
//                 N <= NV+2+NK4*7+(NUM OF KVAL(I)<>0) 
//       IFLAG...INDICATES IF SUCCESSFUL OR NOT, 
//                        = 1 SUCCESSFUL RETURN 
//                        <>1 FAILUE, T,RCOEF ARE INVALD. 
// ***WORK* 
//   WK1(N),WK2(MM)......WHERE MM=MAX(IRC*5+105,N*9) 
void blgcs_(const int *ibci, const double *ts, const double *te, int ncd,
	int nv,const int *kval,const double *val,int nk4,const int *idk,const double *rcir,
	int iv,	int irc, double *wk1, double *wk2, int *n, double *t, double *rcoef, int *iflag
);
