/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGisects_HH_
#define _MGisects_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Pvector.h"
#include "mg/isect.h"

//Forward class declaration.
class MGObject;
class MGCCisect_list;
class MGCSisect_list;
class MGSSisect_list;
class MGCFisect_vector;
class MGFFisect;
class MGHHisect;
class MGHHisect_vector;

///MGisects defines a vector of MGisect.
///The vector is implemeted using MGPvector<MGisect>.
///All the methods to handle the vector are available from the MGPvector class,
///and public member m_is_vector. Refer to MGPvector template class.
///MGisects is used to represent an array of intersection lines of 
///two objects.
///The behavior of MGisects is like an auto_ptr. Copy or assignment
///of MGisects means transfer of the ownership of all the included MGisect
///to copied or assigned MGisects and original MGisects does not have the
///ownership any more. Users should be aware of this fact.
///Intersections are obtained from two objects, which are known using
///the member functions object1() and object2().
///****NOTE****
///When two objects' manifold dimension are the same, object1 is this object
///at the invocation of MGObject::intersection(), and object2 is the argument
///object.
///However, their manifold dimension are not the same, object1 is always
///the lower dimension's object and object2 is the higer dimension's object.
class MGCLASS MGisects{

public:

MGPvector<MGisect> m_is_vector;
typedef MGPvector<MGisect>::iterator              iterator;
typedef MGPvector<MGisect>::const_iterator        const_iterator;
typedef MGPvector<MGisect>::reverse_iterator       reverse_iterator;
typedef MGPvector<MGisect>::const_reverse_iterator const_reverse_iterator;
typedef MGPvector<MGisect>::reference             reference;
typedef MGPvector<MGisect>::const_reference       const_reference;
typedef MGPvector<MGisect>::size_type             size_type;

///String stream Function
MGDECL friend std::ostream& operator << (std::ostream& ostrm, const MGisects& is);

//////////// Constructor ////////////

///Constructor(of size 0)
MGisects(
	const MGObject* obj1=0,
	const MGObject* obj2=0
);

///Construct from MGCCisect_list.
MGisects(MGCCisect_list& ccis);

///Construct from MGCSisect_list.
MGisects(MGCSisect_list& csis);

///Construct from MGCSisect_list.
MGisects(MGSSisect_list& ssis);

///Construct from MGCFisect_vector.
MGisects(MGCFisect_vector& cfis);

///Construct from MGHHisect.
MGisects(MGHHisect& hhi);

///Construct from MGCFisect_vector.
MGisects(MGHHisect_vector& hhis);

//Copy Constructor.
//MGisects(const MGisects& vector);

//////////// Destructor. ////////////
//~MGisects(){;};

//////////// Operator overload. ////////////

///Assignment.
///MGisects& MGisects::operator= (const MGisects&);

const MGisect* operator[](size_t i)const{return m_is_vector[i];};
MGisect* operator[](size_t i){return m_is_vector[i];};

///Assignment.
///MGisects& MGisects::operator= (const MGisects&);

//////////// Member Function. ////////////

///Return(but does not remove) last element in the vector.
///If vector is empty, behavior is undefined.
const MGisect* back() const{return m_is_vector.back();};
MGisect* back() {return m_is_vector.back();};

///Get the iterator of the first element of the m_is_vector.
iterator begin(){return m_is_vector.begin();}
const_iterator begin() const{return m_is_vector.begin();}

///Clear all the elements in m_is_vector.
void clear(){m_is_vector.clear();}

///Return true if there are no items in the vector, false(0) otherwise.
bool empty() const{return m_is_vector.empty();};

///Get the iterator of the next of the last element of the m_is_vector.
iterator end(){return m_is_vector.end();}
const_iterator end() const{return m_is_vector.end();}

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
iterator erase(iterator i){return m_is_vector.erase(i);}

///Exchange first and second order of MGisect.
void exchange12();

/// Return(but does not remove) first element in the vector.
/// If vector is empty, behavior is undefined.
const MGisect* front() const{return m_is_vector.front();};
MGisect* front(){return m_is_vector.front();};

///Insert MGisect at the index position i.
///Transfers the ownership of the isect to this vector.
void insertAt(iterator i, MGisect* isect){
	m_is_vector.insert(i, isect);
};

///Get the 1st object pointer of the i-th intersection.
///Generally objects are different for each intersection.
///Ex. in the case of Shell to Shell intersection, different Face pointer
///will be returned.
const MGObject* object1(size_t i)const;

///Get the 2nd object pointer of the i-th intersection.
///Generally objects are different for each intersection.
///Ex. in the case of Shell to Shell intersection, different Face pointer
///will be returned.
const MGObject* object2(size_t i)const;

///Erase the last element of m_is_vector if not null.
void pop_back(){m_is_vector.pop_back();}

///Adds one MGisect* to the end of the vector.
///Transfers the ownership of the isect to this vector.
void push_back(MGisect* isect){m_is_vector.push_back(isect);};

///append all the member of isects to the end of the vector.
///Transfers the ownership of the isect in isects to this vector.
void push_back(MGisects& isects);

/// Return the number of items that are in the vector.
size_t size() const{return m_is_vector.size();};

private:
	const MGObject* m_object1;///< Object 1.
	const MGObject* m_object2;///< Object 2.

};

/** @} */ // end of IsectContainer group
#endif
