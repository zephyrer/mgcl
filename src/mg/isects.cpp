/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/isects.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/FSurface.h"
#include "topo/CFisect_vector.h"
#include "topo/HHisect_vector.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "topo/FFisect.h"

//MGisects defines a vector of MGisect.
//The vector is implemeted using MGPvector.
//All the methods to handle the vector are available from the MGPvector class,
//and public member m_is_vector. Refer to MGPvector template class.
//MGisects is used to represent an array of intersection lines of 
//two objects.
//The behavior of MGisects is like an auto_ptr. Copy or assignment
//of MGisects means transfer of the ownership of all the included MGisect
//to copied or assigned MGisects and original MGisects does not have the
//ownership any more. Users should be aware of this fact.
//

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////// Constructor //////////

//Void constructor(of size 0)
MGisects::MGisects(
	const MGObject* obj1,
	const MGObject* obj2
):m_object1(obj1), m_object2(obj2){
	if(obj1 && obj2)
		assert(obj1->manifold_dimension()<=obj2->manifold_dimension());
}

//Construct from MGCCisect_list.
MGisects::MGisects(MGCCisect_list& ccis)
:m_object1(ccis.curve1()), m_object2(ccis.curve2()){
	MGCCisect_list::iterator i=ccis.begin(), ie=ccis.end();
	for(; i!=ie; i++){
		m_is_vector.push_back(new MGCCisect(*i));
	}
}

//Construct from MGCSisect_list.
MGisects::MGisects(MGCSisect_list& csis)
:m_object1(csis.curve()), m_object2(0){
	const MGFSurface* srf=csis.surface();
	if(srf)
		m_object2=srf->object_pointer();
	MGCSisect_list::iterator i=csis.begin(), ie=csis.end();
	for(; i!=ie; i++){
		m_is_vector.push_back(new MGCSisect(*i));
	}
}

//Construct from MGSSisect_list.
MGisects::MGisects(MGSSisect_list& ssis)
:m_object1(0), m_object2(0){
	const MGFSurface* srf1=ssis.surface1();
	if(srf1)
		m_object1=srf1->object_pointer();
	const MGFSurface* srf2=ssis.surface2();
	if(srf2)
		m_object2=srf2->object_pointer();
	MGSSisect_list::iterator i=ssis.begin(), ie=ssis.end();
	for(; i!=ie; i++){
		m_is_vector.push_back(new MGSSisect(*i));
	}
}

//Construct from MGCFisect_vector.
MGisects::MGisects(MGCFisect_vector& cfis)
:m_object1(cfis.curve()), m_object2(&(cfis[0].face())){
	MGCFisect_vector::iterator i=cfis.begin(), ie=cfis.end();
	for(; i!=ie; i++){
		m_is_vector.push_back(new MGCFisect(*i));
	}
}

//Construct from MGHHisect.
MGisects::MGisects(MGHHisect& hhi){
	size_t n=hhi.num_of_uvline();
	MGCurve* iline;
	MGFPline uvline1;
	MGFPline uvline2;
	for(size_t i=0; i<n; i++){
		hhi.release_front(iline,uvline1,uvline2);
		m_is_vector.push_back(new MGFFisect(iline,uvline1,uvline2));
	}
}

//Construct from MGHHisect_vector.
MGisects::MGisects(MGHHisect_vector& hhis){
	MGHHisect_vector::iterator i=hhis.begin(), ie=hhis.end();
	for(; i!=ie; i++){
		MGisects isects(*i);
		size_t n=isects.size();
		for(size_t j=0; j<n; j++) m_is_vector.push_back(isects.m_is_vector.release(j));
	}
}

//Debug Function
std::ostream& operator << (std::ostream& ostrm, const MGisects& is){
	ostrm<<"MGisects::"<<&is<<", m_object1="<<is.m_object1<<", m_object2="<<is.m_object2;
	size_t n=is.size();
	ostrm<<", number of isect="<<n<<std::endl;
	MGisects::const_iterator itr; size_t i=0;
	for(itr=is.begin(); itr!=is.end(); itr++) ostrm<<i++<<":"<<(*itr);
	ostrm<<std::endl;
	return ostrm;
}

////////// Member Function. //////////

//Replace first and second order of MGisect.
void MGisects::exchange12(){
	bool change=false;
	if(m_object2){
		int m1=m_object1->manifold_dimension(), m2=m_object2->manifold_dimension();
		if(m1==m2){
			const MGObject* objsave=m_object1;
			m_object1=m_object2;
			m_object2=objsave;
			change=true;
		}
	}
	if(size()&&change){
		iterator i=begin(), ie=end();
		for(; i!=ie; i++) (**i).exchange12();
	}
}

//Get the 1st object pointer of the i-th intersection.
const MGObject* MGisects::object1(size_t i)const{
	return m_is_vector[i]->object1(m_object1);
}

//Get the 2nd object pointer of the i-th intersection.
const MGObject* MGisects::object2(size_t i)const{
	return m_is_vector[i]->object1(m_object2);
}

//append all the member of isects to the end of the vector.
//Transfers the ownership of all the isect in isects to this vector.
void MGisects::push_back(MGisects& isects){
	for(size_t j=0,n=isects.size(); j<n; j++) m_is_vector.push_back(isects.m_is_vector.release(j));
}
