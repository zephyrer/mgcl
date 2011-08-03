/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSSisect_list_HH_
#define _MGSSisect_list_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/SSisect.h"
#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<MGSSisect>;
#pragma warning( pop )
#else
#include <list>
#endif

//Forward class declaration.
class MGFSurface;

/// MGSSisect_list defines linked list of MGSSisect.
/// This is value based list.
/// Used to represent intersection lines of two surfaces.
///The behavior of MGSSisect is like a auto_ptr. Copy or assignment
///of MGSSisect means transfer of the ownership of all the included curve
///to copied or assigned MGSSisect and original MGSSisect does not have the
///curve any more. User should be aware of this fact.
class MGCLASS MGSSisect_list{

public:
#if defined(MGCL_DLL)
	typedef MGListProxy<MGSSisect> container_type;
#else
	typedef std::list<MGSSisect> container_type;
#endif

container_type m_SSilist;

typedef container_type::iterator SSiterator;
typedef container_type::const_iterator const_SSiterator;

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

///String stream Function
MGDECL friend std::ostream & operator << (std::ostream&, const MGSSisect_list& );

//////////// Constructor ////////////
explicit MGSSisect_list(const MGFSurface *s1=NULL, const MGFSurface *s2=NULL)
: m_surface1(s1), m_surface2(s2){;};

///Copy Constructor.
///MGSSisect_list(const MGSSisect_list& list);

//////////// Destructor. ////////////
~MGSSisect_list(){;};

//////////// Operator overload. ////////////

///Assignment.
///MGSSisect_list& MGSSisect_list::operator= (const MGSSisect_list&);

//////////// Member Function. ////////////

/// Adds the MGSSisect to the end of the list.
///isect transfer the ownership of the curves in isect to this list.
void append(const MGSSisect& isect);
void push_back(const MGSSisect& isect){m_SSilist.push_back(isect);};
void append(const MGSSisect_list& isectlist);

/// 全てのコンポーネントを指定して交線を追加
///Add one intersection line to the list.
///iline, param1, and param2 must be newed objects, and their ownership
///are transfered to MGSSisect_list.
void append(
	MGCurve* iline,
	MGCurve* param1,
	MGCurve* param2,
	const MGSSRELATION r1=MGSSREL_UNKNOWN);

/// 全てのコンポーネントを指定して交線を追加
///Add one intersection line to the list.
///this append copies the three curves.
void append(
	const MGCurve& iline,
	const MGCurve& param1,
	const MGCurve& param2,
	const MGSSRELATION r1=MGSSREL_UNKNOWN);

///Get the pointer of the first element of the m_SSilist.
SSiterator begin(){return m_SSilist.begin();}
const_SSiterator begin() const{return m_SSilist.begin();}

///Clear all the elements in m_SSilist.
void clear(){m_SSilist.clear();}

///Get the pointer of the next of the last element of the m_SSilist.
SSiterator end(){return m_SSilist.end();}
const_SSiterator end() const{return m_SSilist.end();}

///Find where in this ssi2  have common parts (in line_zero()) in 
///their world representation.
///Fucntion's return value is the iterator of this that had the common.
///		!=end():have common part. 
///		==end():no common part(except a point) found.
SSiterator find_common(const MGSSisect& ssi2);

///Return the pointer to surface1.
const MGFSurface* surface1() const {return m_surface1;}

///Return the pointer to surface2.
const MGFSurface* surface2() const {return m_surface2;}

/// Return the number of items that are in the list.
size_t entries() const{return m_SSilist.size();};
size_t size() const{return m_SSilist.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
SSiterator erase(SSiterator i){return m_SSilist.erase(i);}

/// Return(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
const MGSSisect& first() const{return m_SSilist.front();};
const MGSSisect& front() const{return m_SSilist.front();};
MGSSisect& front(){return m_SSilist.front();};

///Insert MGSSisect at the index position i.
///This position must be between zero and the number of items in the list,
/// or behavior is undefined.
void insertAt(SSiterator i, const MGSSisect& isect)
{m_SSilist.insert(i, isect);};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
bool isEmpty() const{return m_SSilist.empty();};
bool empty() const{return m_SSilist.empty();};

/// Return(but does not remove) lst element in the list.
/// If list is empty, behavior is undefined.
const MGSSisect& last() const{return m_SSilist.back();};
const MGSSisect& back() const{return m_SSilist.back();};
MGSSisect& back(){return m_SSilist.back();};

///Erase the last element of m_SSilist if not null.
void pop_back(){m_SSilist.pop_back();}

///Erase the first element of m_SSilist if not null.
void pop_front(){m_SSilist.pop_front();}

/// Adds the MGSSisect to the beginning of the list.
///isect transfer the ownership of the curves in isect to this list.
void prepend(const MGSSisect& isect){push_front(isect);};
void push_front(const MGSSisect& isect){m_SSilist.push_front(isect);};

///Remove the MGSSisect and return the MGSSisect. If i is no valid, 
/// behavior is undefined.
MGSSisect removeAt(SSiterator i);

///Remove the first MGSSisect int the list and return the MGSSisect.
///If i is not valid, behavior is undefined.
MGSSisect removeFirst();

///Remove the last MGSSisect in the list and return the MGSSisect.
///If i is not valid, behavior is undefined.
MGSSisect removeLast();

///Replace first and second order of surface 1 and 2.
MGSSisect_list& replace12() ;

private:
	const MGFSurface *m_surface1;	///< Surface 1.
	const MGFSurface *m_surface2;	///< Surface 2.

};

/** @} */ // end of IsectContainer group
#endif
