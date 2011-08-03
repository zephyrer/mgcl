/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/GelPosition.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGGelPosition Class.
//MGGelPosition is a class which expresses which group a gel belongs to.
//

//Assignment.
MGGelPosition& MGGelPosition::operator= (const MGGelPosition& GelPosition2){
	m_gel=GelPosition2.m_gel;m_group=GelPosition2.m_group;
	return *this;
}

//Equal operator
bool MGGelPosition::operator==(const MGGelPosition& gelp2) const{
	if(m_gel!=gelp2.m_gel)
		return false;
	if(m_group!=gelp2.m_group)
		return false;
	return true;
}

bool MGGelPosition::operator<(const MGGelPosition& gp2)const{
	const MGGroup* g1=group();
	const MGGroup* g2=gp2.group();
	if(g1!=g2)
		return g1<g2;

	return m_gel<gp2.m_gel;
}

//String output function.
std::ostream& operator<<(std::ostream& out, const MGGelPosition& pos){
	out<<"MGGelPosition::m_gel="<<pos.m_gel<<",m_group="<<pos.m_group;
	return out;
}

//Return if the bottom gel is a MGGroup or MGObject.
bool MGGelPosition::gel_is_object()const{
	if(m_gel){
		return dynamic_cast<MGObject*>(m_gel)!=0;
	}else
		return false;
}

//Generate a newed clone object.
MGGelPosition* MGGelPosition::clone()const{
	return new MGGelPosition(*this);
}

//perform add operation of this gel position.
//(insert the gel of gelp in the group of gelp)
void MGGelPosition::do_add(){
	group()->push_back(gel());
}

//perform remove operation of this gel position.
//(Release the gel of m_gel from the group of m_group, but does not delete the gel).
void MGGelPosition::do_remove(){
	assert(m_group);
	MGGroup::reverse_iterator i=m_group->rbegin(), iend=m_group->rend(), inext;
	while(i!=iend){
		inext=i; inext++;
		if(*i==m_gel){
			m_group->release(inext.base());//release the found gel
			break;
		}
		i=inext;
	}
}

//Get the object pointer of this.
const MGObject* MGGelPosition::object()const{
	if(m_gel)
		return m_gel->includes_object();
	return 0;
}
MGObject* MGGelPosition::object(){
	if(m_gel)
		return m_gel->includes_object();
	return 0;
}

//Test if this is symmetric to gel2.
//Symmetric means:
//both gels are MGObject and they have the same manifold dimension.
bool MGGelPosition::symmetric(const MGGelPosition& gp2)const{
	const MGGel* gel1=gel();
	const MGObject* obj1=gel1->object();
	if(!obj1)
		return false;
	const MGGel* gel2=gp2.gel();
	const MGObject* obj2=gel2->object();
	if(!obj2)
		return false;
	if(obj1->manifold_dimension()!=obj2->manifold_dimension())
		return false;

	return true;
}
