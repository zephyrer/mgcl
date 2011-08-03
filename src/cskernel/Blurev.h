/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//  BLUREV REVERSES THE DIRECTION OF PARAMETER OF A GIVEN B-REP. 
// *** INPUT  * 
//      K,N1,T1(N1+K),RCOEF1(IRC1,NCD),IRC1,NCD....ORIGINAL B-REP. 
//      IRC2........ROW DIMENSION OF RCOEF2 
// *** OUTPUT * 
//      N2,T2(N2+K),RCOEF2(IRC2,NCD)...NEW B-REP OBTAINED. N2=N1 
// *** NOTE * 
//      N2,T2,RCOEF2 MAY BE THE SAME AREA AS N1,T1,RCOEF1, EACH. 
void blurev_(int k, int n1, const double *t1, 
	const double *rcoef1, int irc1, int ncd, int irc2, 
	int *n2, double *t2, double *rcoef2);

