/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/CCisect_list.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCCisect_list defines linked list of MGCCisect.
// Used to represent intersection points of two curves.

// Constructor
MGCCisect_list::MGCCisect_list(const MGCurve* c1, const MGCurve*c2)
: m_curve1(c1), m_curve2(c2)
{
	double er1=0., er2=0.;
	if(c1) er1=c1->param_error();
	if(c2) er1=c2->param_error();
	m_error=(er1*er1+er2*er2)*9.;	//set error.
}
//Copy Constructor.

// Destructor.
//MGCCisect_list::~MGCCisect_list(){m_curve1=m_curve2=0;}

// Operator overload.

//Assignment.

// Member Function.

void MGCCisect_list::append(const MGCCisect& isect){
// Adds the MGCCisect to the end of the list.
	CCiterator itr; double dif1,dif2;
	for(itr=begin(); itr!=end(); itr++){
		dif1=(*itr).param1()-isect.param1();
		dif2=(*itr).param2()-isect.param2();
		if((dif1*dif1+dif2*dif2)<=m_error) return;
	}

	push_back(isect);
}

// 交点の全てのコンポーネントを指定して，交点リストに追加
//Add one intersection point to the list.
void MGCCisect_list::append(
		const MGPosition& point,	//Intesection point(x,y,)
		double t1,					//parameter value of curve 1.
		double t2,					//parameter value of curve 2.
		const MGCCRELATION r1){
	append(MGCCisect(point,t1,t2,r1));
}

void MGCCisect_list::append(const MGCCisect_list& list){
// Adds the MGCCisect_list to the end of the list.
	const_CCiterator i;
	for(i=list.begin(); i!=list.end(); i++) append(*i);
}

MGCCisect MGCCisect_list::removeAt(CCiterator i){
//Remove the MGCCisect and return the MGCCisect. If i is no valid, 
// behavior is undefined.
	MGCCisect isect=*i;
	erase(i);
	return isect;
}

MGCCisect MGCCisect_list::removeFirst(){
//Remove the first MGCCisect int the list and return the MGCCisect.
//If i is not valid, behavior is undefined.
	MGCCisect isect=front();
	pop_front();
	return isect;
}

MGCCisect MGCCisect_list::removeLast(){
//Remove the first MGCCisect int the list and return the MGCCisect.
//If i is not valid, behavior is undefined.
	MGCCisect isect=back();
	pop_back();
	return isect;
}

MGCCisect_list& MGCCisect_list::replace12() {
//Replace first and second order of curve 1 and 2.

	//Replace curve pointer.
	const MGCurve* save=m_curve1;
	m_curve1=m_curve2; m_curve2=save;

	//Replace parameter 1 and 2.
	double t1,t2;
	CCiterator itr;
	for(itr=begin(); itr!=end(); itr++){
		t1=(*itr).param1(); t2=(*itr).param2();
		(*itr)=MGCCisect((*itr).point(),t2,t1,(*itr).rel());
	}
	return *this;
}
