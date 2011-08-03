/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLPoints_HH_
#define _mgTLPoints_HH_

#include "mg/MGCL.h"
#include <vector>
#include "mg/Position.h"

#if defined(MGCL_DLL)
#pragma warning(push)
#pragma warning(disable : 4231)
#pragma warning(disable : 4660)
MGTEMPLATE template class MGCLASS std::vector<MGPosition>;
#pragma warning(pop)
#endif

///mgTLPoints holds the vector of the surface parameter (u,v).
///ポリゴンの頂点保持クラス.
class MGCLASS mgTLPoints{

public:
	typedef std::vector<MGPosition>::iterator iterator;
	typedef std::vector<MGPosition>::const_iterator const_iterator;
	std::vector<MGPosition> m_positions;	///<ポリゴン頂点のベクトル

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLPoints& tlPoints);

//////////// constructor ///////////////
mgTLPoints():m_positions(){;};

////////// member function /////////////
iterator begin(){return m_positions.begin();}
const_iterator begin()const{return m_positions.begin();}
const_iterator end()const{return m_positions.end();}
iterator end(){return m_positions.end();}

size_t add(double u, double v);
size_t add(const MGPosition &position){ return add(position[0],position[1]);};
size_t size()const{return m_positions.size();};

/*
///Retrieve the normal data of i-th element.
void normal(size_t i, MGVector& N);

///Retrieve (u,v) data of i-th element.
void uv(size_t i, double& u, double& v);
void uv(size_t i, MGPosition& UV);

///Retrieve the world coordinate data of i-th element.
void world(size_t i, MGPosition& P);
*/

MGPosition& operator[](size_t pos){return m_positions[pos];};
const MGPosition& operator[](size_t pos) const{return m_positions[pos];};

private:

};

#endif
