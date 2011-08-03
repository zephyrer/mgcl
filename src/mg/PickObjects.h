/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#pragma once 

#ifndef _MGPickObjects_HH_
#define _MGPickObjects_HH_

#include <vector>
#include "mg/AbstractGels.h"
#include "mg/Pvector.h"
#include "mg/PickObject.h"

class MGCurve;
class MGFSurface;
class MGOfstream;
class MGIfstream;
class MGGelPositions;
class MGPickObject;
class COpenGLWindow;
class CPoint;

/////////////////////////////////////////////////

/** @addtogroup MGObjectRelated
 *  @{
 */

///a container class for MGPickObject.
class MGPickObjects{

public:
	typedef MGPvector<MGPickObject> container_type;

	/// types:
	typedef container_type::reference              reference;
	typedef container_type::const_reference        const_reference;
	typedef container_type::iterator               iterator;
	typedef container_type::const_iterator         const_iterator;
	typedef container_type::size_type              size_type;
	typedef container_type::reverse_iterator       reverse_iterator;
	typedef container_type::const_reverse_iterator const_reverse_iterator;

/// Constructors

MGPickObjects(){;};

///Copy constructor
MGPickObjects(const MGPickObjects& pobjs);

///Construct MGPickObjects of one pobj.
MGPickObjects(
	const MGPickObject& pobj
);

///virtual ~MGPickObjects();

///Operator overload.
MGPickObjects& operator=(const MGPickObjects& pobjs);

///Set operation.
MGPickObjects& operator+=(const MGPickObjects& gelps){push_back(gelps);return *this;};
MGPickObjects& operator+=(const MGPickObject& gelp);
MGPickObjects& operator-=(const MGPickObjects& gelps){remove(gelps);return *this;};
MGPickObjects& operator-=(const MGPickObject& gelp){remove(gelp);return *this;};
MGPickObjects& operator-=(const MGAbstractGels& types){remove(types);return *this;};
MGPickObjects& operator&=(const MGPickObjects& gelps){reset_with_common(gelps);return *this;};

///append the current objects(MGGelPositions).
void append_object(const MGGelPositions& gelps);

///Replace this sequence with [first,last).
void assign(const_iterator first, const_iterator last);

///convert this pick objects to gels. If a face that is a part of a shell,
///the face pointer will be converted to the shell pointer.
void convert_to_ShellGels(MGGelPositions& gels)const;

const MGPickObject& back()const{return *m_PickObjects.back();};
MGPickObject& back(){return *m_PickObjects.back();};
iterator begin(){return m_PickObjects.begin();};
reverse_iterator rbegin(){return m_PickObjects.rbegin();};
const_iterator begin()const{return m_PickObjects.begin();};
const_reverse_iterator rbegin()const{return m_PickObjects.rbegin();};
void clear(){m_PickObjects.clear();};
bool empty()const{return m_PickObjects.empty();};
iterator end(){return m_PickObjects.end();};
reverse_iterator rend(){return m_PickObjects.rend();};
const_iterator end()const{return m_PickObjects.end();};
const_reverse_iterator rend()const{return m_PickObjects.rend();};
const MGPickObject& front()const{return *m_PickObjects.front();};
MGPickObject& front(){return *m_PickObjects.front();};
const MGPickObject& operator[](size_t i)const{return *m_PickObjects[i];};
MGPickObject& operator[](size_t i){return *m_PickObjects[i];};
void pop_back(){m_PickObjects.pop_back();};

///find the same pobj in this objects.
iterator find(const MGPickObject& pobj);
const_iterator find(const MGPickObject& pobj)const;

///Test if there is a MGPickObject whose gel is the input gelin
///in this MGPickObjects' member. If the gel of MGPickObject is MGShell and gelin
///is MGFace, test is performed if the shell includes the face gelin
///as the shell constituent.
///Returns true if gelin is included in this MGPickObjects.
bool includes(const MGGel* gelin)const;

///erase element MGPickObject gelp if found in this.
///Function's return value is the following iterator
///of the erased element x.
iterator erase(const MGPickObject& pobj);

/// erase sequence [first, last).
void erase(iterator first, iterator last);

/// erase sequence i.
iterator erase(iterator i);

/// erase i-th element.
void erase(size_t i){erase(begin()+i);};

/// erase after the elments after the front().
///Resutl has length 1 sequence.
void erase_except_front();

///Make the display list of this object as a highlighted one.
void make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density=1///<line density to draw surface in wire mode.
)const;

///add one pobj.
///Function's return value is the numbe of PickObjects defined.
size_t push_back(const MGPickObject& pobj);
size_t push_back(const MGPickObjects& pobjs);

///Remove pobj if found in this.
void remove(const MGPickObject& pobj){erase(pobj);};
void remove(const MGPickObjects& pobjs);

///Remove objects of type from this pickobjects.
void remove(const MGAbstractGels& types);

///Remove gelps from this pickobjects.
void remove(const MGGelPositions& gelps);

///Select objects of specified type from this and reset with them.
void reset_objects(const MGAbstractGels& types);

///replace this with the common objects of this and pobjs2.
void reset_with_common(const MGPickObjects& pobjs2);

///replace this with symmetric_differecne of this and pobj, that is;
///(1) remove the same MGPickObject from this and pobjs2.
///(2) append the result pobjs2 to this.
void reset_with_symmetric_difference(const MGPickObjects& pobjs2);

///reserve the size n, which are all null.
void reserve(size_t n);

///resize the length of the sequence.
void resize(size_t n){m_PickObjects.resize(n);};

///resize the length of the sequence.
void reset(size_t i, const MGPickObject& pobj);

///Select objects of input type from this.
///Function's return value is pickobjects selected.
///This will be unchanged.
MGPickObjects select(const MGAbstractGels& types)const;

///Select the 1st MGCurve from this.
///Function's return value is MGPickObject of MGCurve 1st encountered in this
///MGPickObject sequence. If this did not includes any MGCurve,
///null MGPickOjbect will be returned.
///This will be unchanged.
MGPickObject select_1st_curve()const;

///Select all the MGCurve from this.
///MGPickObject of MGCurve encountered in this MGPickObject sequence will be appended
///in curves.
///This will be unchanged.
void select_curves(MGPickObjects& curves)const;

///Select the 1st MGFSurface from this.
///Function's return value is MGPickObject of MGFSurface 1st encountered in this
///MGPickObject sequence. If this did not includes any MGFSurface,
///null MGPickObject will be returned.
///This will be unchanged.
MGPickObject select_1st_fsurface()const;

///Select all the MGFSurface from this.
///MGPickObjects of MGFSurface encountered in this MGPickObject sequence will be appended
///in surfaces.
///This will be unchanged.
void select_fsurfaces(MGPickObjects& surfaces)const;

///Obtain the pobj number defined.
size_t size()const{return m_PickObjects.size();};

container_type& object_vector(){return m_PickObjects;};
const container_type& object_vector()const{return m_PickObjects;};

protected:
	container_type m_PickObjects;///vector of MGPickObject.

};

/** @} */ // end of MGObjectRelated group
#endif // _MGPickObjects_HH_
