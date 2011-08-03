/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/FOuterCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGFOuterCurve Class.
//MGFOuterCurve is to represent Face's outer boundary in combination of
//loops and perimeter curves. This is private class for Face and Loop.

///////Operator oveload///////

///////Member function///////

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGFOuterCurve& focrv){
	if(focrv.is_loop()){
		out<<"m_loop="<<focrv.m_loop;
	}else{
		out<<"m_peri_id="<<focrv.m_peri_id<<", m_t=("
			<<focrv.m_t[0]<<","<<focrv.m_t[1]<<")";
	}
	return out;
}
