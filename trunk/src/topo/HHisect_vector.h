/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGHHisect_vector_HH_
#define _MGHHisect_vector_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include <vector>
#include "topo/HHisect.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGHHisect>;
#pragma warning( pop )
#endif

//Forward class declaration.

///MGHHisect_vector defines a vector of MGHHisect.
///The vector is implemeted using STL's vector.
///All the methods to handle the vector are available from the STL's vector class,
///and public member m_HHivector. Refer to STL vector class.
///MGHHisect_vector is used to represent intersection lines of a shell with
///another shell, a face, or a surface.
///The behavior of MGHHisect is like an auto_ptr. Copy or assignment
///of MGHHisect means transfer of the ownership of all the included curves
///to copied or assigned MGHHisect and original MGHHisect does not have the
///ownership more. Users should be aware of this fact.
class MGCLASS MGHHisect_vector{

public:

std::vector<MGHHisect> m_HHivector;

typedef std::vector<MGHHisect>::iterator HHiterator;
typedef std::vector<MGHHisect>::const_iterator const_HHiterator;

typedef std::vector<MGHHisect>::iterator iterator;
typedef std::vector<MGHHisect>::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream & operator << (std::ostream&, const MGHHisect_vector& );

//////////// Constructor ////////////

///Void constructor.
MGHHisect_vector(){;};

///Constructor of 1 MGHHisect.
MGHHisect_vector(MGHHisect& hhi):m_HHivector(1,hhi){;};

///Copy Constructor.
///MGHHisect_vector(const MGHHisect_vector& vector);

//////////// Destructor. ////////////
///~MGHHisect_vector(){;};

//////////// Operator overload. ////////////

///Assignment.
///MGHHisect_vector& MGHHisect_vector::operator= (const MGHHisect_vector&);

const MGHHisect& operator[](size_t i)const{return m_HHivector[i];};
MGHHisect& operator[](size_t i){return m_HHivector[i];};

//////////// Member Function. ////////////

///Adds one MGHHisect to the end of the vector.
///Transfers the ownership of the curves in isect to this vector.
void push_back(MGHHisect& isect){m_HHivector.push_back(isect);};
void push_back(MGHHisect_vector& isects);

///Return(but does not remove) last element in the vector.
///If vector is empty, behavior is undefined.
const MGHHisect& back() const{return m_HHivector.back();};
MGHHisect& back() {return m_HHivector.back();};

///Get the pointer of the first element of the m_HHivector.
HHiterator begin(){return m_HHivector.begin();}
const_HHiterator begin() const{return m_HHivector.begin();}

///Clear all the elements in m_HHivector.
void clear(){m_HHivector.clear();}

///Get the pointer of the next of the last element of the m_HHivector.
HHiterator end(){return m_HHivector.end();}
const_HHiterator end() const{return m_HHivector.end();}

/// Return the number of items that are in the vector.
size_t size() const{return m_HHivector.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
HHiterator erase(HHiterator i){return m_HHivector.erase(i);}

/// Return(but does not remove) first element in the vector.
/// If vector is empty, behavior is undefined.
const MGHHisect& front() const{return m_HHivector.front();};
MGHHisect& front(){return m_HHivector.front();};

///Insert MGHHisect at the index position i.
void insertAt(HHiterator i, MGHHisect& isect){
	m_HHivector.insert(i, isect);
};

///Return true if there are no items in the vector, false(0) otherwise.
bool empty() const{return m_HHivector.empty();};

///Erase the last element of m_HHivector if not null.
void pop_back(){m_HHivector.pop_back();}

///Replace first and second order of surface 1 and 2.
MGHHisect_vector& replace12();

};

/** @} */ // end of IsectContainer group
#endif
