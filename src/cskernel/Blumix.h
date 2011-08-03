/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUMIX WILL MIX TWO TWO DIMENSIONAL B-REP INTO ONE THREE D B-REP. 
// *** INPUT * 
//  ERROR......ERROR ESTIMATE OF BLIP 
//  KCOD1,K1,N1,T1(N1+K1),RCOEF(IRC1,2).....DESCRIBE THE 1ST B-REP 
//           OF ORDER K1. KCOD1 IS COORDINATE KIND OF THE B-REP, MUST BE 
//           1,2, OR,3 (YZ,ZX,XY EACH) 
//  KCOD2,K2,N2,T2(N2+K2),RCOE2(IRC2,2).....DESCRIBE THE SECOND B-REP OF 
//           ORDER K2. KCOD2 MAY BE 4, INDICATING GIRTH B-REP, I.E. 
//           RCOE2(.,1) IS THE DATA POINTS OF THE 1ST B-REP 
//           AND RCOE2(.,2) IS THE MISSING COORDINATE OF THE 1ST B-REP. 
//  IRC......ROW DIMENSION OF THE VARIABLE RCOEF, I.E. LENGTH OF RCOEF 
//           ALLOWED TO USE. 
//           IRC MUST BE GREATER THAN 11 AND N1. 
// *** OUTPUT * 
//  VECS(3),VECE(3).......START AND END 1ST DERIVATIVES. 
//  NOUT,T(N),RCOEF(IRC,3)......THREE DIMENSIONAL DATA POINTS 
//                           (ORDINATE AND ABSCISSA). 
//  IFLAG........OUTPUT IF SUCCESSFUL RETURN OR NOT. 
//           =0: SUCCESSFUL RETURN 
//           =2: MEMORY EXAUSTED. 
// ***WORK AREA*** 
//  KSEQ(IRC),WORK(NN)....WHERE NN=MAX(IRC*4,4*K*K+3*K).(K=MAX(K1,K2)) 
void blumix_(double error, int kcod1, int k1, 
	int n1,const double *t1,const double *rcoe1, int irc1, 
	int kcod2, int k2, int n2,const double *t2,const double *rcoe2,
	int irc2, int irc, int *kseq, double *work, 
	double *vecs, double *vece, int *nout, double *t, 
	double *rcoef, int *iflag);
