/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/LSPoint.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGLSPoint Class.
//MGLSPoint is to express a loop and a surface intersection point.
//The expression is {MGEdge* binder, double tb,(u,v)}, where binder is
//binder edge of the loop, tb is parameter value of the binder, and (u,v) is
//the surface parameter value.

///////Constructor////////

///////Operator oveload///////

//Comparison operator.
bool MGLSPoint::operator< (const MGLSPoint& ls2)const{
	if(m_pedge!=ls2.m_pedge) return (*m_pedge)<(*ls2.m_pedge);
	return m_t<ls2.m_t;
}

bool MGLSPoint::operator== (const MGLSPoint& ls2)const{
	if(m_pedge!=ls2.m_pedge) return false;
	if(!MGREqual(m_t, ls2.m_t)) return false;
	if(!MGREqual(m_u, ls2.m_u)) return false;
	return MGREqual(m_v, ls2.m_v);
}

///////Member function///////

void MGLSPoint::set_surface_param(const MGPosition& uv){
	m_u=uv[0]; m_v=uv[1];
}

//Return surface's parameter data.
MGPosition MGLSPoint::surface_param()const{return MGPosition(m_u, m_v);}
void MGLSPoint::surface_param(double& u, double& v)const{u=m_u; v=m_v;}

//Obtain world point coordinate data from the binder edge.
MGVector MGLSPoint::world_point()const{
	const MGEdge* bedge=parameter_edge()->binder_edge();
	return bedge->eval(m_t);
}

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGLSPoint& lsp){
	out<<"MGLSPoint::m_pedge="<<lsp.m_pedge<<", m_t="<<lsp.m_t;
	out<<"(m_u, m_v)=("<<lsp.m_u<<","<<lsp.m_v<<")";
	return out;
}
