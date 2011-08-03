/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BVUTAN GETS APPROXIMATE UNIT TANGENT VECTOR, GIVEN POINT SEQUENCE. 
// THE VECTOR MAY BE STARTING AT P(1,.), OR ENDING AT P(NP,.) 
// ***INPUT* 
//   ISE....INDICATES WHETHER STARTING OR ENDING VECTOR 
//           =1: STARTING AT P(1,.)    =2: ENDING AT P(NP,.) 
//   NP....GIVES NUM OF INPUT POINTS 
//   P(IP,NCD)....PROVIDES POINT SEQUENCE OF NCD SPACE DIMENSION 
//   IP....IS ROW DIMENSION OF THE VARIABLE P 
//   NCD...IS SPACE DIMENSION OF TNPUT POINTS 
// ***OUTUT* 
//   RTAN(NCD).....THE UNIT TANGENT VECTOR 
//   IFLAG....=1: SUCCESSFUL RETURN 
//            =2: SOME ERROR DETECTED AND RTAN(.) WAS NOT CALCULATED 
// ***WORK* 
//   WORK(81+NCD*8) 
void bvutan_(int ise, int np,const double *p, 
	int ip, int ncd, double *work, double *rtan, 
	int *iflag);