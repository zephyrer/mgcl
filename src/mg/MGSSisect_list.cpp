/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/FSurface.h"
#include "mg/SSisect_list.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSSisect_list defines singly linked list of MGSSisect.
// Used to represent intersection points of two surfaces.

//	list<MGSSisect> ; // List of intersection points. 
//	MGSurface *m_surface1;	// Surface 1.
//	MGSurface *m_surface2;	// Surface 2.

// Constructor

//Copy Constructor.

// Destructor.

// Operator overload.

//Assignment.

// Member Function.

void MGSSisect_list::append(const MGSSisect& isect){
// Adds the MGSSisect to the end of the list.
	push_back(isect);
}
void MGSSisect_list::append(const MGSSisect_list& isectlist){
	const_SSiterator i=isectlist.begin(), ie=isectlist.end();
	for(; i!=ie; i++)
		append(*i);
}

// 全てのコンポーネントを指定して交線を追加
//Add one intersection line to the list.
//iline, param1, and param2 must be newed objects, and their ownership
//are transfered to MGSSisect_list.
void MGSSisect_list::append(
	MGCurve* iline,
	MGCurve* param1,
	MGCurve* param2,
	const MGSSRELATION r1){
// Adds the MGSSisect to the end of the list.
	push_back(MGSSisect(iline,param1,param2,r1));
}

// 全てのコンポーネントを指定して交線を追加
//Add one intersection line to the list.
void MGSSisect_list::append(
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION r1){
// Adds the MGSSisect to the end of the list.
	push_back(MGSSisect(iline,param1,param2,r1));
}

//Find where in this ssi2  have common parts (in line_zero()) in 
//their world representation.
//Fucntion's return value is the iterator of this that had the common.
//		!=end():have common part. 
//		==end():no common part(except a point) found.
MGSSisect_list::SSiterator MGSSisect_list::find_common(const MGSSisect& ssi2){
	SSiterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if(ssi2.has_common(*i)) return i;
	}
	return ie;
}

MGSSisect MGSSisect_list::removeAt(SSiterator i){
//Remove the MGSSisect and return the MGSSisect. If i is no valid, 
// behavior is undefined.
	MGSSisect isect=(*i);
	erase(i);
	return isect;
}

MGSSisect MGSSisect_list::removeFirst(){
//Remove the first MGSSisect int the list and return the MGSSisect.
//If i is not valid, behavior is undefined.
	MGSSisect isect=front();
	pop_front();
	return isect;
}

MGSSisect MGSSisect_list::removeLast(){
//Remove the first MGSSisect int the list and return the MGSSisect.
//If i is not valid, behavior is undefined.
	MGSSisect isect=back();
	pop_back();
	return isect;
}

MGSSisect_list& MGSSisect_list::replace12() {
//Replace first and second order of surface 1 and 2.

	//Replace surface pointer.
	const MGFSurface* save=m_surface1;
	m_surface1=m_surface2; m_surface2=save;

	//Replace intersection param line.
	SSiterator i=begin(), ie=end();
	for(; i!=ie; i++) i->replace12();
	return *this;
}
