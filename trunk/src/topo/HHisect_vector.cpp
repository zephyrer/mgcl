/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/HHisect_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGHHisect_vector defines a vector of MGHHisect.
//The vector is implemeted using STL's vector.
//All the methods to handle the vector are available from the STL's vector class,
//and public member m_HHivector. Refer to STL vector class.
//MGHHisect_vector is used to represent intersection lines of a shell with
//another shell, a face, or a surface.
//The behavior of MGHHisect is like an auto_ptr. Copy or assignment
//of MGHHisect means transfer of the ownership of all the included curves
//to copied or assigned MGHHisect and original MGHHisect does not have the
//ownership more. Users should be aware of this fact.
//

////////// Constructor //////////

////////// Operator overload. //////////

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGHHisect_vector& hhis){
	size_t n=hhis.size();
	out<<"MGHHisect_vector::number of isect="<<n<<std::endl;
	MGHHisect_vector::const_HHiterator i;
	size_t j=0;
	for(i=hhis.begin(); i!=hhis.end(); i++) out<<j++<<":"<<(*i);
	return out;
}

//Assignment.
//MGHHisect_vector& MGHHisect_vector::operator= (const MGHHisect_vector&);

////////// Member Function. //////////
void MGHHisect_vector::push_back(MGHHisect_vector& isects){
	size_t n=isects.size();
	for(size_t i=0; i<n; i++){
		push_back(isects[i]);
	}
	isects.clear();
}

//Replace first and second order of surface 1 and 2.
MGHHisect_vector& MGHHisect_vector::replace12(){
	MGHHisect_vector::HHiterator i=begin(), ie=end();
	for(; i!=ie; i++) (*i).exchange12();
	return *this;
}
