/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/GelPositions.h"
#include "mg/PickObjects.h"
#include "mg/Gel.h"
#include "mg/isects.h"
#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "topo/Face.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGGelPositions Class.

//////////Constructor//////////

//Copy constructor.
MGGelPositions::MGGelPositions(const MGGelPositions& gpos2)
:m_gelps(gpos2.m_gelps){
}

//Conversion constructor of MGPickObjects.
MGGelPositions::MGGelPositions(const MGPickObjects& gelp){
	size_t n=gelp.size();
	for(size_t i=0; i<n; i++){
		MGGelPosition geli=gelp[i];
		m_gelps.push_back(geli);
	}
}

///////////////operator overloaded//////////////

//////////Member Function//////////

// Output virtual function.
std::ostream& operator<<(std::ostream& ostrm, const MGGelPositions& gelps){
	ostrm<<"MGGelPositions::number of gelps="<<gelps.size()<<std::endl;
	MGGelPositions::const_iterator i=gelps.begin(), ie=gelps.end();	
	for(size_t j=0; i!=ie; i++, j++){
		ostrm<<"gelp-"<<j<<":"<<(*i)<<std::endl;
	}
	return ostrm;
}

//Replace this sequence with [first,last).
void MGGelPositions::assign(const_iterator first, const_iterator last){
	MGGelPositions pobjs;
	for(;first!=last; first++)
		pobjs.push_back(*first);
	m_gelps.clear();
	push_back(pobjs);
}

// erase element MGGelPosition gelp. Function's return value is the following iterator
// of the erased element x.
MGGelPositions::iterator MGGelPositions::erase(const MGGelPosition& gelp){
	iterator i=find(gelp), ie=end();
	if(i!=end()) return erase(i);
	return ie;
}

//Test if there is a MGGelPosition whose gel is the input gelin
//in this GelPositions' member. If the gel of MGGelPosition is MGShell and gelin
//is MGFace, test is performed if the shell includes the face gelin
//as the shell constituent.
//Returns true if gelin is included in this MGGelPositions.
bool MGGelPositions::includes(const MGGel* gelin)const{
	size_t ngelp=size();
	for(size_t i=0; i<ngelp; i++){
		const MGGel* geli=m_gelps[i].gel();
		if(geli == gelin)
			return true;

		const MGShell* shl=geli->shell();
		const MGFace* f=gelin->face();
		if(!shl || !f)
			continue;

		if(shl->includes(f))
			return true;
	}
	return false;
}

//Remove objects of type from this.
void MGGelPositions::remove(const MGAbstractGels& types){
	int n=size();
	for(int i=n-1; i>=0; i--){
		if(m_gelps[i].gel()->type_is(types))
			m_gelps.erase(m_gelps.begin()+i);
	}
}

//Remove gelps from this.
void MGGelPositions::remove(const MGGelPositions& gelps){
	size_t ngelp=gelps.size();
	for(size_t i=0; i<ngelp; i++){
		const MGGelPosition& po=gelps[i];
		iterator hit = find(po);
		if(hit != end()){//if po was already in current object, erase it.
			erase(hit);	
		}
	}
}
//reserve the size n, which are all null.
void MGGelPositions::reserve(size_t n){
	m_gelps.reserve(n);
}

void MGGelPositions::reset(size_t i, const MGGelPosition& pobj){
	m_gelps[i]=pobj;
}

//Select objects of specified type from this and reset with them.
void MGGelPositions::reset_objects(const MGAbstractGels& types){
	int n=size();
	for(int i=n-1; i>=0; i--){
		if(!((*this)[i].gel()->type_is(types))) erase(begin()+i);
	}
}

//replace this with the common objects of this and pobjs2.
void MGGelPositions::reset_with_common(const MGGelPositions& pobjs2){
	const_iterator je=pobjs2.end();
	for(int i=size()-1; i>=0; i--){
		if(pobjs2.find((*this)[i])==je)
			erase(i);
	}
}

//replace this with symmetric_differecne of this and pobj, that is;
//(1) remove the same MGGelPosition from this and pobjss.
//(2) append the result pobjs2 to this.
//On return, pobjs2 will have null sequence.
void MGGelPositions::reset_with_symmetric_difference(MGGelPositions& pobjs2){
	int n1=size(), n2=pobjs2.size();
	MGGelPositions* p1=this;
	MGGelPositions* p2=&pobjs2;
	int n=n1;
	if(n1>n2){
		n=n2;
		p1=&pobjs2; p2=this;
	}
	for(int i=n-1; i>=0; i--){
		iterator j=p2->find((*p1)[i]);
		if(j!=p2->end()){
			p2->erase(j);
			p1->erase(p1->begin()+i);
		}
	}
	push_back(pobjs2);
}

//Select objects of input type from this.
//Function's return value is MGGelPositions selected.
//This will be unchanged.
MGGelPositions MGGelPositions::select(const MGAbstractGels& types)const{
	MGGelPositions pobjs2;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((*i).gel()->type_is(types))
			pobjs2.push_back(*i);
	}
	return pobjs2;
}

//Select the 1st MGCurve from this.
//Function's return value is MGGelPosition of MGCurve 1st encountered in this
//MGGelPosition sequence. If this did not includes any MGCurve,
//null MGGelPosition will be returned.
//This will be unchanged.
MGGelPosition MGGelPositions::select_1st_curve()const{
	const MGCurve* curve=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGGel* obj=(*i).object();
		curve=dynamic_cast<const MGCurve*>(obj);
		if(curve)
			return *i;
	}
	return MGGelPosition();
}

//Select all the MGCurve from this.
//MGGelPosition of MGCurve encountered in this MGGelPosition sequence will be appended
//in curves.
//This will be unchanged.
void MGGelPositions::select_curves(MGGelPositions& curves)const{
	const MGCurve* curve=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGGel* obj=(*i).object();
		curve=dynamic_cast<const MGCurve*>(obj);
		if(curve){
			curves.push_back(*i);
		}
	}
}

//Select the 1st MGFSurface from this.
//Function's return value is MGGelPosition of MGFSurface 1st encountered in this
//MGGelPosition sequence. If this did not includes any MGFSurface,
//null MGGelPosition will be returned.
//This will be unchanged.
MGGelPosition MGGelPositions::select_1st_fsurface()const{
	const MGFSurface* f=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGGel* obj=(*i).object();
		f=dynamic_cast<const MGFSurface*>(obj);
		if(f)
			return *i;
	}
	return MGGelPosition();
}

//Select all the MGFSurface from this.
//MGGelPositions of MGFSurface encountered in this MGGelPosition sequence will be appended
//in surfaces.
//This will be unchanged.
void MGGelPositions::select_fsurfaces(MGGelPositions& surfaces)const{
	const MGFSurface* surface;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGGel* obj=(*i).object();
		surface=dynamic_cast<const MGFSurface*>(obj);
		if(surface){
			surfaces.push_back(*i);
		}
	}
}

//Test if this is symmetric to gels2.
//Symmetric means:
//(1) number of gels included is the same.
//(2) all of the gels are MGObject and they have the same manifold dimension.
bool MGGelPositions::symmetric(
	 const MGGelPositions& gels2
)const{
	size_t n=size();
	if(gels2.size()!=n)
		return false;

	for(size_t i=0; i<n; i++){
		if(m_gelps[i].symmetric(gels2[i]))
			continue;
		return false;
	}
	return true;
}
