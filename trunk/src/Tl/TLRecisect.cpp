/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "Tl/TLRecisect.h"
#include "Tl/TLRect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//mgTLRecisect is a proprietry class for Face tessellation.
//mgTLRecisect represents one intesection of Loop and a subdivided (u,v)
//rectangle's perimeter.

///////Operator oveload///////
bool mgTLRecisect::operator< (const mgTLRecisect& is2)const{
	bool in1=going_in(), in2=is2.going_in();
	if(in1 && !in2) return true;
	if(!in1 && in2) return false;
	if(m_perim<is2.m_perim) return true;
	if(m_perim>is2.m_perim) return false;
	if(m_perim<=1) return (*m_is).t()<(*(is2.m_is)).t();
	return (*m_is).t()>(*(is2.m_is)).t();
}

bool mgTLRecisect::operator== (const mgTLRecisect& is)const{
	if(m_perim!=is.m_perim) return false;
	return (*m_is).t()==(*(is.m_is)).t();
}

bool mgTLRecisect::going_in()const{
	if(m_perim==0 || m_perim==3) return m_is->increase()==1;
	return m_is->increase()==-1;
}

bool mgTLRecisect::going_out()const{
	if(m_perim==0 || m_perim==3) return m_is->increase()==-1;
	return m_is->increase()==1;
}

//Test if this recisect is in lower range or upper range when mgTLRect is
//subdivided at s.
bool mgTLRecisect::greater_than(
	bool is_u,		//input if s is u value(true), or not.
	double s		//parameter value on a perimeter.
){
	if(is_u){//If s is u value.
		if(m_perim==1) return true;
		if(m_perim==3) return false;
	}else{
		if(m_perim==2) return true;
		if(m_perim==0) return false;
	}
	return (*m_is).t()>s;
}

//Test if ri2's intersection loop and this are the same.
bool mgTLRecisect::same_loop(const mgTLRecisect& ri2) const{
	const MGLoop* lp1=m_is->isect().loop();
	const MGLoop* lp2=ri2.m_is->isect().loop();
	return lp1==lp2;
}

void mgTLRecisect::set_going_in(){
	if(m_perim==0 || m_perim==3) m_is->set_increase(1);
	else m_is->set_increase(-1);
}
void mgTLRecisect::set_going_out(){
	if(m_perim==0 || m_perim==3) m_is->set_increase(-1);
	else m_is->set_increase(1);
}

//Set the flag that this isect is processed by createTrimPolygon().
void mgTLRecisect::set_processed(){
	m_perim=-(m_perim+1);
}

//Return the status of this rec isects.
// =0:unknown; =1:going in, =2:going out,
int mgTLRecisect::status()const{
	if(going_in()) return 1;
	if(going_out()) return 2;
	return 0;
}

//Get the (u,v) parameter value of this intersection point.
MGPosition mgTLRecisect::uv(
	const mgTLRect& rect
)const{
	switch(perim()){
	case 0: return MGPosition(is()->t(), rect.vmin());
	case 1: return MGPosition(rect.umax(), is()->t());
	case 2: return MGPosition(is()->t(), rect.vmax());
	case 3: return MGPosition(rect.umin(), is()->t());
	default: return MGPosition();
	}
}

ostream& operator<< (ostream& out, const mgTLRecisect& ris){
	out<<"TLRecisect::m_perim="<<ris.m_perim;
	if(ris.m_perim!=-9999)
		out<<",m_is=("<<*(ris.m_is)<<")";
	return out;
}
