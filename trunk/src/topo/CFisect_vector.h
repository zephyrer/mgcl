/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCFisect_vector_HH_
#define _MGCFisect_vector_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/CFisect.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGCFisect>;
#pragma warning( pop )
#endif

///MGCFisect_vector defines a vector of MGCFisect.
///The vector is implemeted using STL's vector.
///All the methods to handle the vector are available from the STL's vector class,
///and public member m_CFivector. Refer to STL vector class.
///MGCFisect_vector is used to represent intersection points of a curve with a shell.
///The behavior of MGCFisect is like an auto_ptr. Copy or assignment
///of MGCFisect means transfer of the ownership of all the included curves
///to copied or assigned MGCFisect and original MGCFisect does not have the
///ownership more. Users should be aware of this fact.
class MGCLASS MGCFisect_vector{

public:

std::vector<MGCFisect> m_CFivector;

typedef std::vector<MGCFisect>::iterator iterator;
typedef std::vector<MGCFisect>::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream & operator << (std::ostream&, const MGCFisect_vector& );

//////////// Constructor ////////////

///Void constructor.
MGCFisect_vector(const MGCurve* curve=0):m_curve(curve){;};

///Constructor of 1 MGCFisect.
MGCFisect_vector(const MGCurve* curve,MGCFisect& cfi):m_curve(curve),m_CFivector(1,cfi){;};

///Copy Constructor.
///MGCFisect_vector(const MGCFisect_vector& vector);

//////////// Destructor. ////////////
///~MGCFisect_vector(){;};

//////////// Operator overload. ////////////

///Assignment.
///MGCFisect_vector& MGCFisect_vector::operator= (const MGCFisect_vector&);

const MGCFisect& operator[](size_t i)const{return m_CFivector[i];};
MGCFisect& operator[](size_t i){return m_CFivector[i];};

//////////// Member Function. ////////////

///Adds one MGCFisect to the end of the vector.
///Transfers the ownership of the curves in isect to this vector.
void push_back(MGCFisect& isect){m_CFivector.push_back(isect);};

///Return(but does not remove) last element in the vector.
///If vector is empty, behavior is undefined.
const MGCFisect& back() const{return m_CFivector.back();};
MGCFisect& back() {return m_CFivector.back();};

///Get the pointer of the first element of the m_CFivector.
iterator begin(){return m_CFivector.begin();}
const_iterator begin() const{return m_CFivector.begin();}

///Clear all the elements in m_CFivector.
void clear(){m_CFivector.clear();}

///Get the curve pointer.
const MGCurve* curve()const{return m_curve;};

///Get the pointer of the next of the last element of the m_CFivector.
iterator end(){return m_CFivector.end();}
const_iterator end() const{return m_CFivector.end();}

/// Return the number of items that are in the vector.
size_t size() const{return m_CFivector.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
iterator erase(iterator i){return m_CFivector.erase(i);}

/// Return(but does not remove) first element in the vector.
/// If vector is empty, behavior is undefined.
const MGCFisect& front() const{return m_CFivector.front();};
MGCFisect& front(){return m_CFivector.front();};

///Insert MGCFisect at the index position i.
void insertAt(iterator i, MGCFisect& isect){
	m_CFivector.insert(i, isect);
};

///Return true if there are no items in the vector, false(0) otherwise.
bool empty()const{return m_CFivector.empty();};

///Erase the last element of m_CFivector if not null.
void pop_back(){m_CFivector.pop_back();}

private:
	const MGCurve* m_curve;///<Curve pointer of the curve to get the intersection.

};

/** @} */ // end of IsectContainer group
#endif
