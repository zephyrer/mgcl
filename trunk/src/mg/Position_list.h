/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPosition_list_HH_
#define _MGPosition_list_HH_
/** @addtogroup IsectContainer
 *  @{
 */

#include "mg/Position.h"
#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<MGPosition>;
#pragma warning( pop )
#else
#include <list>
#endif

//Forward class declaration.
class MGBox;
class MGCurve;
class MGFSurface;
class MGIfstream;
class MGOfstream;
class MGSSisect;

///MGPosition_list provides a list of Positions.
///MGPosition_list is used to store curve or surface parameter as
///two objects intersection data.
class MGCLASS MGPosition_list{

public:
#if defined(MGCL_DLL)
	typedef MGListProxy<MGPosition> container_type;
#else
	typedef std::list<MGPosition> container_type;
#endif

container_type m_Plist;

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;

//////////// Constructor ////////////

/// Void constructor(Constructor of length 0).
MGPosition_list(){;};

/// Constructor of length 1.
MGPosition_list(const MGPosition& P): m_Plist(1,P){;};

///Construct MGPosition_list by replacing each element P of list by Q,
///where, Q=MGPosition(P.sdim(),P,start1,start2).
MGPosition_list(const MGPosition_list& list,	///Original MGPosition_list 
		size_t start1,	///Start position to store of elements Q of
						///new MGPosition_list.
		size_t start2);	///Start position to retrieve of elements P of list.

///Copy Constructor.
///MGPosition_list(const MGPosition_list& list);

////////////// Destructor. ////////////
~MGPosition_list(){;};

//////////// Operator overload. ////////////

///Debug Function
MGDECL friend std::ostream& operator << (std::ostream&, const MGPosition_list& );

///Assignment.
///MGPosition_list& MGPosition_list::operator= (const MGPosition_list&);

//////////// Member Function. ////////////

/// Adds the parameter to the end of the list.
///Function's return value is true if appended, false if not.
bool append(const MGPosition& pos);
void push_back(const MGPosition& pos){m_Plist.push_back(pos);};

///Add parameter pair p of crv1 and 2 to the end of list if p is not
///included in the list already.
///Function's return value is true if appended, false if not.
bool append(
	const MGCurve& crv1,
	const MGCurve& crv2,
	const MGPosition& p
);

///Add parameter all the data in the list to the end of list if the member p is not
///included in the list already.
void append(
	const MGCurve& crv1,
	const MGCurve& crv2,
	const MGPosition_list& list
);

///Add parameter uv of surface to the end of list if uv is not
///included in the list already.
///Function's return value is true if appended, false if not.
bool append(
	const MGFSurface& srf,
	const MGPosition& uv
);

///Add parameter all the uv's of list to the end of this list if the
///uv is not included in this list already.
void append(
	const MGFSurface& srf,
	const MGPosition_list& list
);

///Add parameter pair tuv of crv and surface to the end of list if tuv is not
///included in the list already.
bool append(
	const MGCurve& crv,
	const MGFSurface& srf,
	const MGPosition& tuv
);

///Add parameter pair uvuv of surface srf1 and sur2 to the end of list
/// if uvuv is not included in the list already.
///Function's return value is true if appended, false if not.
bool append(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
){ return add(true,srf1,srf2,uvuv);}

///Add parameter pair uvuv of surface srf1 and sur2 to the end of list
void append(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition_list& list
);

///Get the pointer of the first element of the m_Plist.
iterator begin(){return m_Plist.begin();}
const_iterator begin() const{return m_Plist.begin();}

///Clear all the elements in m_Plist.
void clear(){m_Plist.clear();}

///Get the pointer of the next of the last element of the m_Plist.
iterator end(){return m_Plist.end();}
const_iterator end() const{return m_Plist.end();}

/// Returns the number of items that are in the list.
size_t entries() const{	return m_Plist.size();};
size_t size() const{return m_Plist.size();};

///Erase the element of iterator i.
///Returned is the iterator located after the element i.
iterator erase(iterator i){return m_Plist.erase(i);}

/// Returns(but does not remove) first element in the list.
/// If list is empty, behavior is undefined.
const MGPosition& first() const{return m_Plist.front();};
const MGPosition& front() const{return m_Plist.front();};
MGPosition& first(){return m_Plist.front();};
MGPosition& front(){return m_Plist.front();};

///Test if one of the points in the list is included in
///the box. If so, return id of the point. The id is of the first point
///found in the list.
///If no points was not included in the list, end() wil be returned in id.
///Function's return value is true if a point is in the box, and false
///if no points are in the box.
///n is the number of space dimension of points in this list to check.
///If n>0 is specified, the first n coordinates of the position data in
///the list are checked if included in the box.
///If n==0 is specified, all the coordinate data are checked.
bool in(const MGBox& box, const_iterator& id, size_t n=0) const;
bool in(const MGBox& box, iterator& id, size_t n=0);

///Inserts position at the index i.
///This index must be between zero and the number of items in the list,
/// or behavior is undefined.
void insertAt(iterator i, const MGPosition& pos)
{m_Plist.insert(i, pos);};

///Return true (1) if there are no items in the list,
/// false(0) otherwise.
bool isEmpty() const{return m_Plist.empty();};
bool empty() const{return m_Plist.empty();};

/// Return(but does not remove) last element in the list.
/// If list is empty, behavior is undefined.
const MGPosition& last() const{return m_Plist.back();};
const MGPosition& back() const{return m_Plist.back();};
MGPosition& last(){return m_Plist.back();};
MGPosition& back(){return m_Plist.back();};

///Erase the last element of m_Plist if not null.
void pop_back(){m_Plist.pop_back();}

///Erase the first element of m_Plist if not null.
void pop_front(){m_Plist.pop_front();}

/// Adds the parameter to the beginning of the list.
void prepend(const MGPosition& pos){m_Plist.push_front(pos);};
void push_front(const MGPosition& pos){m_Plist.push_front(pos);};

///Add parameter pair uvuv of surface srf1 and sur2 to the beginning of list
///if uvuv is not included in the list already.
bool prepend(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
){ return add(false,srf1,srf2,uvuv);};

///Remove position in the list that is the same point as P.
///When distace of the two point is within error, they are regarded as same.
///Function's return value is num of points deleted.
int remove(
	double error,		///<square of error allowed to regard same.
	const MGPosition& P,///<Position to remove.
	size_t n=0			///<Number of space dimension to check,
						///<If n==0, all the coordinates are checked.
);

///Remove parameter pair uv of surface srf1 from the list
/// if uv is included in the list.
int remove(
	const MGFSurface& srf,
	const MGPosition& uv
);

///Remove parameter pair uvuv of surface srf1 and sur2 from the list
/// if uvuv is included in the list.
///uvuv's space dimension is 4. And uvuv(0) and uvuv(1) are srf1's
/// parameter (u,v), uvuv(2) and uvuv(3) are srf2's parameter (u,v).
///Function's return value is true if removal was done, and false if
///no same point was found and removal was not performed.
int remove(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
);

///Remove same elements in this list(parameter pair uvuv of surface srf1 and sur2)
///if elements are within the tolerance.
///Space dimension of the elements(uvuv) is 4. And uvuv(0) and uvuv(1) are srf1's
/// parameter (u,v), uvuv(2) and uvuv(3) are srf2's parameter (u,v).
///Function's return value is number of removed elements.
int remove_uvuv(
	const MGFSurface& srf1,
	const MGFSurface& srf2
);

///Remove the position and return the position. If i is not valid, 
/// behavior is undefined.
MGPosition removeAt(iterator i);

///Remove the first position in the list and return the position.
///If i is not valid, behavior is undefined.
MGPosition removeFirst();

///Remove the last position in the list and return the position.
///If i is not valid, behavior is undefined.
MGPosition removeLast();

///Remove uvuv that is on the ssi.
///Function's return value is the number of points deleted.
size_t removeOn(
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGSSisect& ssi
);

///reverse the ordering of the elements in the list.
void reverse_order();

///Sort the positions in the list according to the surface parameter space
///ordering. Positions (uvuv(id),uvuv(id+1)) in the list is a parameter
///position of a face.
void sort_uv_space(size_t id);

private:

///Dump Functions
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );

///Add parameter pair uvuv of surface srf1 and sur2 to the end(append==true)
///or the beginning of list if uvuv is not included in the list already.
bool add(
	bool append,			
	const MGFSurface& srf1,
	const MGFSurface& srf2,
	const MGPosition& uvuv
);

};

/** @} */ // end of IsectContainer group
#endif
