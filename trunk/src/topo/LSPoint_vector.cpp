/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Surface.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/LSPoint_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLSPoint_vector defines a vector of MGLSPoint.
// Used to represent intersection points of a loop and a surface.

//////////Constructor///////////

//Copy Constructor.
//MGLSPoint_vector(const MGLSPoint_vector& vec);

// Destructor.
//~MGLSPoint_vector(){;};

// Operator overload.

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGLSPoint_vector& vec){
	out<<"MGLSPoint_vector::";
	size_t n=vec.entries();
	out<<"number of points="<<n<<std::endl;
	MGLSPoint_vector::const_LSiterator itr; size_t i=0;
	for(itr=vec.m_lspoints.begin(); itr!=vec.m_lspoints.end(); itr++)
		out<<i++<<":"<<(*itr);
	return out;
}

//Assignment.
//MGLSPoint_vector& MGLSPoint_vector::operator= (const MGLSPoint_vector&);

// Member Function.

//Add one intersection point to the list.
bool MGLSPoint_vector::append(const MGLSPoint& lsp){
// Adds the MGLSPoint to the end of the list.
	double errorSqr=MGTolerance::line_zero()*1.5; errorSqr*=errorSqr;
	MGVector P=lsp.world_point();
	LSiterator itr=m_lspoints.begin(), itrend=m_lspoints.end();
	for(; itr!=itrend; itr++){
		if(lsp.parameter_edge()!=itr->parameter_edge()) continue;
		MGVector Q=P-(itr->world_point());
		if(Q%Q<=errorSqr) return false;
	}
	m_lspoints.push_back(lsp);
	return true;
}

// Adds the MGLSPoint_vector to the end of the vector.
void MGLSPoint_vector::append(const MGLSPoint_vector& vec){
// Adds the MGLLisect_vector to the end of the list.
	const_LSiterator i;
	for(i=vec.m_lspoints.begin(); i!=vec.m_lspoints.end(); i++) append(*i);
}
