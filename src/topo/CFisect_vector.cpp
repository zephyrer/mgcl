/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/CFisect_vector.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGCFisect_vector defines a vector of MGCFisect.
//The vector is implemeted using STL's vector.
//All the methods to handle the vector are available from the STL's vector class,
//and public member m_CFivector. Refer to STL vector class.
//MGCFisect_vector is used to represent intersection lines of a shell with
//another shell, a face, or a surface.
//The behavior of MGCFisect is like an auto_ptr. Copy or assignment
//of MGCFisect means transfer of the ownership of all the included curves
//to copied or assigned MGCFisect and original MGCFisect does not have the
//ownership more. Users should be aware of this fact.
//

////////// Constructor //////////

////////// Operator overload. //////////

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGCFisect_vector& cfis){
	size_t n=cfis.size();
	out<<"MGCFisect_vector::number of isect="<<n<<std::endl;
	MGCFisect_vector::const_iterator i;
	size_t j=0;
	for(i=cfis.begin(); i!=cfis.end(); i++) out<<j++<<":"<<(*i);
	return out;
}

//Assignment.
//MGCFisect_vector& MGCFisect_vector::operator= (const MGCFisect_vector&);

////////// Member Function. //////////

