/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLTriangles_HH_
#define _mgTLTriangles_HH_

#include "mg/Pvector.h"
#include "Tl/TLTriangle.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<mgTLTriangle*>;
MGTEMPLATE template class MGCLASS MGPvector<mgTLTriangle>;
#pragma warning( pop )
#endif

class MGFace;
class mgTLData;

///A vector of mgTriangle's.
///mgTriangle holds multiple triangles, and mgTLTriangles holds vector of mgTriangle's.
///ポリゴンのベクトル保持クラス
class MGCLASS mgTLTriangles{

public:
	typedef MGPvector<mgTLTriangle>::const_iterator const_iterator;
	typedef MGPvector<mgTLTriangle>::iterator iterator;
	typedef MGPvector<mgTLTriangle>::const_iterator const_triIterator;
	typedef MGPvector<mgTLTriangle>::iterator triIterator;
	MGPvector<mgTLTriangle> m_triangles;	///ポリゴンのベクトル

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLTriangles& tlTriangles);

//////////// constructor ///////////////

mgTLTriangles():m_triangles(){;}

//////////// operator overload ////////////

///Return the i-th mgTLTriangle.
const mgTLTriangle& operator[](size_t i)const;
mgTLTriangle& operator[](size_t i);

////////// member function /////////////
triIterator begin(){return m_triangles.begin();}
triIterator end(){return m_triangles.end();}
const_triIterator begin()const{return m_triangles.begin();}
const_triIterator end()const{return m_triangles.end();}
///const MGFace* getFace()const{return m_f;};

mgTLTriangle* front(){return m_triangles.front();};
void print(std::ostream& out, const mgTLData& tld)const;
void push_back(mgTLTriangle* pTlTriangle){m_triangles.push_back(pTlTriangle);};
void push_back(mgTLTriangles& tris){m_triangles.push_back(tris.m_triangles);};
size_t size()const{return m_triangles.size();};

private:
};

#endif
