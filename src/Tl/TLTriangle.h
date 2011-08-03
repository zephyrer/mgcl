/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLTriangle_HH_
#define _mgTLTriangle_HH_

#include "mg/MGCL.h"
#include <vector>

//#if defined(MGCL_DLL)
//#pragma warning( push )
//#pragma warning( disable : 4231 )
//MGTEMPLATE template class MGCLASS std::vector<size_t>;
//#pragma warning( pop )
//#endif

class MGPosition;
class mgTLPoints;
class MGSurface;
class MGFace;
class mgTLRect;
class mgTLData;

///mgTLTriangle holds (multiple) triangles data, which are indices of mgTLPoints.
///mgTLTriangle is generated from trimmed mgTLRect's using utility class mgTLFan,
///or mgTLFans.
///３角形ポリゴンのタイプ(triFan,triStrip)と頂点インデックスベクトル保持クラス.
class MGCLASS mgTLTriangle{

public:

	typedef std::vector<size_t>::iterator IndexItr;
	typedef std::vector<size_t>::const_iterator CIndexItr;
	typedef std::vector<size_t>::iterator iterator;
	typedef std::vector<size_t>::const_iterator const_iterator;
	std::vector<size_t> m_indices;	///<三角形頂点情報インデックスのベクトル
		///<Let id=m_indices[i],
		///<then mgTLPoints[id]=MGPosition& uv=(u,v) parameter value of the surface.

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLTriangle& triangle);

//////////// constructor ///////////////
mgTLTriangle(
	mgTLRect* rect,	///<mgTLRect this belongs to.
	mgTESTRIANG type = mgTESTRIANG_UNKNOWN		///<３角形頂点タイプのリスト
):m_type(type),m_rect(rect){;};

///construct a triangle whose type is mgTESTRIANG_FAN, and number of
///vertices are 3 including the center.
mgTLTriangle(
	mgTLRect* rect,	///<mgTLRect this belongs to.
	size_t center, size_t v1, size_t v2	///<３角形頂点
);

///construct a triangle whose type is mgTESTRIANG_FAN, and number of
///vertices is n.
mgTLTriangle(
	mgTLRect* rect,	///<mgTLRect this belongs to.
	size_t n,			///<number of vertices.
	const size_t* vertices	///<vertices of array length n.
);

//////////// Operator overload ///////////////
size_t operator[](size_t i)const{return m_indices[i];};

//////////// Member Function ///////////////
CIndexItr begin()const{return m_indices.begin();}
CIndexItr end()const{return m_indices.end();}
IndexItr begin(){return m_indices.begin();}
IndexItr end(){return m_indices.end();}

///Erase i-th element of m_indices.
IndexItr erase(size_t i);

///タイプを返却する mgTESTRIANG_FAN mgTESTRIANG_STRIP
mgTESTRIANG getGeometryType()const{return m_type;};

void push_back(size_t index){m_indices.push_back(index);};
void pop_back(){m_indices.pop_back();};

void print(std::ostream& out, const mgTLData& tld)const;

mgTLRect* rect(){return m_rect;};
const mgTLRect* rect()const{return m_rect;};

///タイプを設定する
void setGeometryType(mgTESTRIANG type){m_type = type;};

///Obtain the number of points included.
size_t size()const{return m_indices.size();};

///Assumed this type is mgTESTRIANG_FAN and at least one of the
///perimeter is textued, compute all the other perimeter's texture.
bool texture_triangle(mgTLData& tldata);

///Compute 3rd texture, given ids of the points.
MGPosition third_texture(
	mgTLData& tldata,
	bool right_hand,///<true if (i1,i2,i) constitues right hand system.
	size_t i1,	///<texture-known point 1.
	size_t i2,	///<texture-known point 2.
	size_t i	///<point to compute texture.
)const;

///Get parameter (u,v) of the surface from the id of this triangle.
const MGPosition& uv(size_t id, const mgTLPoints& tlpoints) const;
MGPosition& uv(size_t id, mgTLPoints& tlpoints);

///Get world coordinates of the surface from the id of this triangl.
///The output is surf.eval(uv(id,tlpoints));
MGPosition world(
	size_t id,	///<id of this triangle. 0<= id <= size().
	const mgTLPoints& tlpoints,///<tlpoints[(*this)[id]] is (u,v),
		///< or uv(id, tlpoints) is (u,v).
	const MGFace& f		///<Tessellated face.
) const;
MGPosition world(
	size_t id,	///<id of this triangle. 0<= id <= size().
	const mgTLPoints& tlpoints,///<tlpoints[(*this)[id]] is (u,v),
		///< or uv(id, tlpoints) is (u,v).
	const MGSurface& surf	///<Tessellated surface. When the object was
		///<a face f, surf=*(f.surface());
) const;

private:
	mgTLRect* m_rect;	///<mgTLRect this triangle belongs to when m_type==mgTESTRIANG_FAN;
				///<when m_type==mgTESTRIANG_STRIP, the 1st mgTLRect neighboring is registered.
	mgTESTRIANG m_type;		///<mgTESTRIANG_FAN, mgTESTRIANG_STRIP

///propagate texture computation from the perimeter (0, textured_vertex)
///to (0,1) perimeter(when increase=false) or to (0,n-1) perimeter
///(when increase=true). The vertex 0 must be already textured.
void propagate_texture(
	mgTLData& tldata,
	int textured_vertex,
	bool increase
);

};

#endif
