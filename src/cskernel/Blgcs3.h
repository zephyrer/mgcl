/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLGCS3         BLGCS3 FOR BLGCS, RESOLVE KNUCKLE INF KVALI(.). 
// *** INPUT * 
//       IBCI(2)....SPECIFY BOUNDARY COND, WHEN KVALI(1)=0,KVALI(NVI)=0. 
//              =1: 1ST DERIV     =2: 2ND DERIV     =3: NO B.C. 
//              (1ST OR 2ND DERIV. PROVIDED IN TS(.) OR TE(.) ) 
//       TS(NCD),TE(NCD)....NECESSARY ONLY WHEN KVALI(.)=0 AND 
//              IBCI(.)=1 OR 2, AND GIVES 1ST OR 2ND DERIV AS ABOVE. 
//       NCD....SPACE DIMENSION OF THE INPUT POINTS, 1<= NCD <=3. 
//       NVI....NUMBER OF INPUT POINTS, I.E. 
//              ( TAUI(I),KVALI(I),VALI(I,.) )  FOR 1<=I<=NVI 
//       TAUI(NVI),KVALI(NVI),VALI(IVI,NCD)....INPUT POINTS WITH KL INF 
//              ( TAUI(I),KVALI(I),VALI(I,.) ) IS ONE PAIR OF THE INPUT. 
//              TAUI(I) IS A DATA POINT AND VALI(I,.), THE CORRESPND DATA. 
//              KVALI(I) IS THE KNUCKLE INF, AND VALI(.,.) IS THE POINT 
//              (I.E. POSITIONAL DATA OF NCD SPACE DIMENSION) 
//       IVI,IV....ROW DIMENSIONS OF THE VARIABLE VALI AND VAL ,EACH. 
// *** OUTPUT * 
//       IBC(2)....GIVE BOUNDARY CONDITIONS FOR BLG4SQ 
//                 (1)=BEGINNING COND, (2)=ENDING COND. 
//                 =1: 1ST DERIV SPECIFIED 
//                 =2: 2ND DERIV SPECIFIED 
//                 =3: NO BOUNDARY COND. 
//       N,TAU(N),VAL(IV,NCD)....GENERATED DATA SEQUENCE INCLUDING 
//              DERIVATIVE INF. THESE CAN BE STRAIGHT INPUT OF BLG4SQ. 
//                N: NUMBER OF DATA 
//                TAU(I),VAL(I,.):  1<=I<=N ARE DATA POINTS AND 
//                    ASSOCIATED DATA. 
//       IFLAG       INDICATES IF SUCCESSFUL OR NOT, 
//                        = 1 SUCCESSFUL RETURN 
//                        <>1 FAILUE, VAL AND TAU ARE INVALID. 
// ***WORK* 
//   WORK(105),KWORK(N) WORK ARRAY OF EACH LENGTH 
//   STATE NUM.(ABOVE L) ARE AS FOLLOWS 
//      0->0,2->0,1->0 STATE =0           1->1,2->1,0->1 STATE =1 
//      0->2 STATE           =2           1->2,2->2 STATE      =3 
//      1->-1 STATE          =4 
//      X->-1                =5 (X<>1) 
void blgcs3_(const int *ibci,const double *ts,const double *te, 
	int ncd, int nvi,const double *taui,const int *kvali, 
	const double *vali, int ivi, int iv, double *work, 
	int *kwork, int *ibc, int *n, double *tau, double *val, int *iflag
);