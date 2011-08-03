/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#include "mg/Curve.h"
#include "mg/FSurface.h"
#include "mg/PickObject.h"
#include "mg/PickObjects.h"
#include "mg/GelPositions.h"
#include "topo/Shell.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////
//a container class for MGPickObject.

//Construct MGPickObjects of one pobj.
MGPickObjects::MGPickObjects(
	const MGPickObject& pobj
):m_PickObjects(pobj.clone()){;}

//Copy constructor
MGPickObjects::MGPickObjects(const MGPickObjects& pobjs){
	size_t n=pobjs.size();
	for(size_t i=0; i<n; i++)
		push_back(pobjs[i]);
}

MGPickObjects& MGPickObjects::operator+=(const MGPickObject& gelp){
	m_PickObjects.push_back(gelp.clone());
	return *this;
}

//Operator overload.
MGPickObjects& MGPickObjects::operator=(const MGPickObjects& pobjs){
	m_PickObjects.clear();
	push_back(pobjs);
	return *this;
}

//append the current objects(MGGelPositions).
void MGPickObjects::append_object(const MGGelPositions& gelps){
	size_t n=gelps.size();
	for(size_t i=0; i<n; i++)
		m_PickObjects.push_back(new MGPickObject(gelps[i]));
}

//Replace this sequence with [first,last).
void MGPickObjects::assign(const_iterator first, const_iterator last){
	MGPvector<MGPickObject> tempobjs;
	for(const_iterator i=first; i!=last; i++)
		tempobjs.push_back((**i).clone());
	m_PickObjects=tempobjs;
}

//convert this pick objects to gels. If a face that is a part of a shell,
//the face pointer will be converted to the shell pointer.
void MGPickObjects::convert_to_ShellGels(MGGelPositions& gels)const{
	gels.clear();
	typedef std::vector<MGShell*> shelVec;
	shelVec shells;
	size_t n=m_PickObjects.size();
	for(size_t i=0; i<n; i++){
		MGPickObject pobji=*(m_PickObjects[i]);
		MGFace* f=pobji.gel()->face();
		if(f){
			MGShell* shl=f->shell();
			if(shl){
				shelVec::iterator found=std::find(shells.begin(), shells.end(),shl);
				if(found!=shells.end())//if found.
					continue;
				pobji.set_gel(shl);
				shells.push_back(shl);
			}
		}
		gels.push_back(pobji);
	}
}

// erase element MGGelPosition gelp. Function's return value is the following iterator
// of the erased element x.
MGPickObjects::iterator MGPickObjects::erase(const MGPickObject& pobj){
	iterator i=find(pobj), ie=end();
	if(i!=end()) return erase(i);
	return ie;
}

// erase sequence [first, last).
void MGPickObjects::erase(iterator first, iterator last){
	m_PickObjects.erase(first,last);
}
MGPickObjects::iterator MGPickObjects::erase(iterator i){
	return m_PickObjects.erase(i);
}

// erase after the elments after the front().
//Resutl has length 1 sequence.
void MGPickObjects::erase_except_front(){
	if(size()<=1) return;
	iterator first=begin()+1, last=end();
	erase(first,last);
}

//find the same pobj in this objects.
MGPickObjects::iterator MGPickObjects::find(
	const MGPickObject& pobj
){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((**i).gel()==pobj.gel())
			return i;
	}
	return ie;
}

//find the same pobj in this objects.
MGPickObjects::const_iterator MGPickObjects::find(
	const MGPickObject& pobj
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((**i).gel()==pobj.gel())
			return i;
	}
	return ie;
}

//Test if there is a MGPickObject whose gel is the input gelin
//in this MGPickObjects' member. If the gel of MGPickObject is MGShell and gelin
//is MGFace, test is performed if the shell includes the face gelin
//as the shell constituent.
//Returns true if gelin is included in this MGPickObjects.
bool MGPickObjects::includes(const MGGel* gelin)const{
	size_t nobj=size();
	for(size_t i=0; i<nobj; i++){
		const MGGel* geli=m_PickObjects[i]->gel();
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

///Make the display list of this object as a highlighted one.
void MGPickObjects::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		(**i).make_display_list_to_hilight(span_length,line_density);
}

//add one pobj.
//Function's return value is the numbe of PickObjects defined.
size_t MGPickObjects::push_back(const MGPickObject& pobj){
	m_PickObjects.push_back(pobj.clone());
	return m_PickObjects.size();
}

size_t MGPickObjects::push_back(const MGPickObjects& pobjs){
	size_t n=pobjs.size();
	for(size_t i=0; i<n; i++)
		push_back(pobjs[i]);
	return m_PickObjects.size();
}

//Remove objects of type from this pickobjects.
void MGPickObjects::remove(const MGAbstractGels& types){
	reverse_iterator i=rbegin(), ie=rend();
	for(; i!=ie;){
		if((**i).gel()->type_is(types)){
			m_PickObjects.erase((++i).base());
		}else
			i++;
	}
}
void MGPickObjects::remove(const MGPickObjects& pobjs){
	reverse_iterator i=rbegin(), iend=rend();
	for(; i!=iend;){
		MGPickObject& po=**i++;
		iterator ifoward=i.base();
		if(pobjs.includes(po.gel()))
			erase(ifoward);	
	}
}

//Remove gelps from this pickobjects.
void MGPickObjects::remove(const MGGelPositions& gelps){
	reverse_iterator i=rbegin(), iend=rend();
	for(; i!=iend;){
		MGPickObject& po=**i++;
		iterator ifoward=i.base();
		if(gelps.includes(po.gel()))
			erase(ifoward);	
	}
}

//reserve the size n, which are all null.
void MGPickObjects::reserve(size_t n){
	m_PickObjects.reserve(n);
}

void MGPickObjects::reset(size_t i, const MGPickObject& pobj){
	m_PickObjects.reset(i,pobj.clone());
}

//Select objects of specified type from this and reset with them.
void MGPickObjects::reset_objects(const MGAbstractGels& types){
	int n=size();
	for(int i=n-1; i>=0; i--){
		if(!(m_PickObjects[i]->gel()->type_is(types))) erase(begin()+i);
	}
}

//replace this with the common objects of this and pobjs2.
void MGPickObjects::reset_with_common(const MGPickObjects& pobjs2){
	const_iterator je=pobjs2.end();
	for(int i=size()-1; i>=0; i--){
		if(pobjs2.find(*m_PickObjects[i])==je)
			erase(i);
	}
}

//replace this with symmetric_differecne of this and pobj, that is;
//(1) remove the same MGPickObject from this and pobjs2.
//(2) append the result pobjs2 to this.
void MGPickObjects::reset_with_symmetric_difference(const MGPickObjects& pobjs2){
	int n1=size(), n2=pobjs2.size();
	MGPickObjects* p1=this;
	MGPickObjects p2contens(pobjs2);
	MGPickObjects* p2=&p2contens;
	int n=n1;
	if(n1>n2){
		n=n2;
		p1=p2; p2=this;
	}
	for(int i=n-1; i>=0; i--){
		iterator j=p2->find((*p1)[i]);
		if(j!=p2->end()){
			p2->erase(j);
			p1->erase(p1->begin()+i);
		}
	}
	push_back(p2contens);
}

//Select objects of input type from this.
//Function's return value is pickobjects selected.
//This will be unchanged.
MGPickObjects MGPickObjects::select(const MGAbstractGels& types)const{
	MGPickObjects pobjs2;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGPickObject* pobji=*i;
		if(pobji->gel()->type_is(types))
			pobjs2.push_back(*pobji);
	}
	return pobjs2;
}

//Select the 1st MGCurve from this.
//Function's return value is MGPickObject of MGCurve 1st encountered in this
//MGPickObject sequence. If this did not includes any MGCurve,
//null MGPickOjbect will be returned.
//This will be unchanged.
MGPickObject MGPickObjects::select_1st_curve()const{
	const MGCurve* curve=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).object();
		curve=dynamic_cast<const MGCurve*>(obj);
		if(curve)
			return **i;
	}
	return MGPickObject();
}

//Select all the MGCurve from this.
//MGPickObject of MGCurve encountered in this MGPickObject sequence will be appended
//in curves.
//This will be unchanged.
void MGPickObjects::select_curves(MGPickObjects& curves)const{
	const MGCurve* curve=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).object();
		curve=dynamic_cast<const MGCurve*>(obj);
		if(curve){
			curves.push_back(**i);
		}
	}
}

//Select the 1st MGFSurface from this.
//Function's return value is MGPickObject of MGFSurface 1st encountered in this
//MGPickObject sequence. If this did not includes any MGFSurface,
//null MGPickObject will be returned.
//This will be unchanged.
MGPickObject MGPickObjects::select_1st_fsurface()const{
	const MGFSurface* f=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).object();
		f=dynamic_cast<const MGFSurface*>(obj);
		if(f)
			return **i;
	}
	return MGPickObject();
}

//Select all the MGFSurface from this.
//MGPickObjects of MGFSurface encountered in this MGPickObject sequence will be appended
//in surfaces.
//This will be unchanged.
void MGPickObjects::select_fsurfaces(MGPickObjects& surfaces)const{
	const MGFSurface* surface;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		const MGObject* obj=(**i).object();
		surface=dynamic_cast<const MGFSurface*>(obj);
		if(surface){
			surfaces.push_back(**i);
		}
	}
}
