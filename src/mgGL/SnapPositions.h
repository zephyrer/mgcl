/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPositions_HH_
#define _MGPositions_HH_

#include <vector>
#include <list>
#include "mg/Position.h"
#include "mg/PickObjects.h"

//
//Define MGSnapPositions Class.

class MGCurve;
class MGGroup;
class MGGel;
class MGFace;
class MGShell;

///MGSnapPositions is a class to store array(vector) of MGPosition's.
class MGCLASS MGSnapPositions{

public:

enum snap_kind{
	DELETE_POINT=-1,
	nopos=0,
	endpos,
	knotpos,
	vertexpos,
	nearpos,
	centerpos,
	ON_SURFACE
};

////////////Constructor////////////

///Void constructor.
///MGSnapPositions();

MGSnapPositions::MGSnapPositions(snap_kind kind=nopos):m_snap_kind(kind){;};

///Virtual Destructor
///virtual ~MGSnapPositions();

////////////Member Function////////////

typedef std::vector<MGPosition>::const_iterator const_iterator;
typedef std::vector<MGPosition>::iterator iterator;
typedef std::pair<const MGObject*, int> obj_num;

///append this positions in m_positions into points.
void append_position_data(std::vector<MGPosition>& points)const;

///Display list name.
size_t dlist_name()const{return size_t(this);};

///Extract position data.
void extract(
	const MGCurve& crv	///<the curve to extract.
);
void extract(
	const MGSurface& srf///<the surface to extract.
);
void extract(
	const MGPoint& point	///<the curve to extract.
);
void MGSnapPositions::extract(
	const MGFace& face	///<the curve to extract.
);
void MGSnapPositions::extract(
	const MGShell& shell	///<the surface to extract.
);
void extract(
	const std::list<MGGel*>& gel_list	///<the list to extract.
);
void extract(
	const MGGel& gel	///<the gel to extract.
);
void extract(
	const MGPickObjects& pobjs	///<array of pick objects to extract.
);

void generate_point_display_list()const;

void get_pick_data(
	const unsigned int* selectBuf, 
	MGPosition& point,	///<point data will be output.
	const MGObject*& obj,///<point's object will be rturned.
	double& t	///<When obj is an MGCurve and this snap kind is nearpos, end, or knot,
		///<the point's parameter value of the curve will be returned.
)const;

///Get and set the snap_kind.
snap_kind get_snap_kind()const{return m_snap_kind;};
void set_snap_kind(snap_kind kind){m_snap_kind=kind;};

const MGPosition& back()const{return m_positions.back();};
MGPosition& back(){return m_positions.back();};
iterator begin(){return m_positions.begin();};
const_iterator begin()const{return m_positions.begin();};
void clear(){m_positions.clear();m_obj_nums.clear();};
bool empty()const{return m_positions.empty();};
iterator end(){return m_positions.end();};
const_iterator end()const{return m_positions.end();};
const MGPosition& front()const{return m_positions.front();};
MGPosition& front(){return m_positions.front();};
const MGPosition& operator[](size_t i)const{return m_positions[i];};
MGPosition& operator[](size_t i){return m_positions[i];};
void pop_back(){m_positions.pop_back();};
void push_back(const MGPosition& pos){m_positions.push_back(pos);};
size_t size(){return m_positions.size();};

///Get the point data array.
const std::vector<MGPosition>& points()const{return m_positions;};
std::vector<MGPosition>& points(){return m_positions;};

private:
std::vector<obj_num> m_obj_nums;///<obj_num includes how many data are included in
	///<m_positions for an MGObject.
	///<m_positions.size()=sum of(m_obj_nums[i].second) for i=0,...,m_obj_nums.size()-1.
std::vector<MGPosition> m_positions;///<All the position data will be stored.
snap_kind m_snap_kind;

};

#endif
