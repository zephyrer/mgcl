/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCCisect_list_HH_
#define _MGCCisect_list_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/CCisect.h"
#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<MGCCisect>;
#pragma warning( pop )
#else
#include <list>
#endif

//Forward class declaration.
class MGCurve;

/// MGCCisect_list defines linked list of MGCCisect.
/// Used to represent Intersection points of a curve 1 and a curve 2.
class MGCLASS MGCCisect_list{

public:

#if defined(MGCL_DLL)
	typedef MGListProxy<MGCCisect> container_type;
#else
	typedef std::list<MGCCisect> container_type;
#endif

container_type m_CCilist;

typedef container_type::iterator CCiterator;
typedef container_type::const_iterator const_CCiterator;

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGCCisect_list& );

/// Constructor
explicit MGCCisect_list(const MGCurve *c1=NULL, const MGCurve *c2=NULL);

//Copy Constructor.
//MGCCisect_list(const MGCCisect_list& list);

//////////// Destructor.//////////
~MGCCisect_list(){;};

////////// Operator overload.//////////

//Assignment.
//MGCCisect_list& MGCCisect_list::operator= (const MGCCisect_list&);

////////////Member Function.////////////

/// Adds the MGCCisect to the end of the list.
void append(const MGCCisect& isect);
void push_back(const MGCCisect& isect){m_CCilist.push_back(isect);};

/// 交点の全てのコンポーネントを指定して，交点リストに追加
///Add one intersection point to the list.
void append(
	const MGPosition& point,	///<Intesection point(x,y,z)
	double t1,					///<parameter value of curve 1.
	double t2,					///<parameter value of curve 2.
	const MGCCRELATION r1=MGCCREL_UNKNOWN
);

/// Adds the MGCCisect_list to the end of the list.
void append(const MGCCisect_list& lst);

///Get the pointer of the first element of the m_CCilist.
CCiterator begin(){return m_CCilist.begin();}
const_CCiterator begin() const{return m_CCilist.begin();}

///Clear all the elements in m_CCilist.
void clear(){m_CCilist.clear();}

///Return the pointer to curve1.
const MGCurve* curve1() const {return m_curve1;}

const MGCurve* curve2() const {return m_curve2;}
///Return the pointer to curve2.

///Get the pointer of the next of the last element of the m_CCilist.
CCiterator end(){return m_CCilist.end();}
const_CCiterator end() const{return m_CCilist.end();}

/// Return the number of items that are in the list.
size_t entries() const{return m_CCilist.size();};
size_t size() const{return m_CCilist.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
CCiterator erase(CCiterator i){return m_CCilist.erase(i);}

/// Return(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
const MGCCisect& first() const{return m_CCilist.front();};
const MGCCisect& front() const{return m_CCilist.front();};

///Insert MGCCisect at the iterator i.
void insertAt(CCiterator i, const MGCCisect& isect)
{m_CCilist.insert(i, isect);};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
bool isEmpty() const{return m_CCilist.empty();};
bool empty() const{return m_CCilist.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGCCisect& last() const{return m_CCilist.back();};
const MGCCisect& back() const{return m_CCilist.back();};

///Erase the first element of m_CCilist if not null.
void pop_front(){m_CCilist.pop_front();}

///Erase the last element of m_CCilist if not null.
void pop_back(){m_CCilist.pop_back();}

/// Adds the MGCCisect to the beginning of the list.
void prepend(const MGCCisect& isect){m_CCilist.push_front(isect);};
void push_front(const MGCCisect& isect){m_CCilist.push_front(isect);};

///Remove the MGCCisect and return the MGCCisect. If i is no valid, 
/// behavior is undefined.
MGCCisect removeAt(CCiterator i);

///Remove the first MGCCisect int the list and return the MGCCisect.
///If i is not valid, behavior is undefined.
MGCCisect removeFirst();

///Remove the last MGCCisect in the list and return the MGCCisect.
///If i is not valid, behavior is undefined.
MGCCisect removeLast();

MGCCisect_list& replace12() ;
///Replace first and second order of curve 1 and 2.

private:
	const MGCurve* m_curve1;	///< Curve 1.
	const MGCurve* m_curve2;	///< Curve 2.
	double m_error;			///< Square of Tolerance in parameter space.

};

/** @} */ // end of IsectContainer group
#endif
