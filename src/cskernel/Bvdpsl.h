/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// FUNCTION TO OBTAIN DISTANCE BETWEEN a LINE(G) AND a POINT(P). 
// INPUT *** 
//    NCD......SPACE DIMENSION OF G AND P, MUST BE 2 OR 3. 
//    G(NCD,2) : PARAMETER OF S.L. AS BELOW 
//            G(.,1):A POINT ON THE S.L., G(.,2):DIRECTIONAL COSINE 
//    P(NCD) : COORDINATE OF A POINT 
// OUTPUT *** BVDPSL : DISTANCE 
double bvdpsl_(int ncd,const double *g,const double *p);
