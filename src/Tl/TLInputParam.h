/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLInputParam_HH_
#define _mgTLInputParam_HH_

////////////

class MGIfstream;
class MGOfstream;
#include "mg/MGCL.h"
#include "mg/Object.h"

///A class that contains all the necessary input parameters to make tessellation.
///This is used to construct mgTLData(the tessellation of a surface), or for
///other parameter for tessellation.
class MGCLASS mgTLInputParam{

public:

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLInputParam& para);

/// Serialization fucntion.
friend MGOfstream& operator<< (MGOfstream& buf, const mgTLInputParam& para);
friend MGIfstream& operator>> (MGIfstream& buf, mgTLInputParam& para);

mgTLInputParam(
	double crvTol=.1,			///<バウンダリのトレランス
	double surfTol=.1,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

///Construct from the object box data and the span length to draw object.
///span_length=MGOpenGLView::span_length().
mgTLInputParam(
	const MGObject& obj,
	double span_length
);

///////////// Operator overload./////////

////////// member function /////////////
double crvTol()const{return m_crvTol;};
void set_crvTol(double new_tol){ m_crvTol = new_tol;}

double surfTol()const{return m_surfTol;};
void set_surfTol(double new_tol){ m_surfTol = new_tol;}

double max_ratio()const{return m_max_ratio;};
void set_max_ratio(double new_ratio){ m_max_ratio = new_ratio;}

MGCL::fan_kind fanKind()const{return m_fk;};
void set_fanKind(MGCL::fan_kind new_fan){ m_fk = new_fan;}

size_t minimum_tri()const{return m_minimum_tri;};
void set_minimum_tri(size_t new_tri){ m_minimum_tri = new_tri;}

double max_edge_len()const{return m_max_edge_len;};
void set_max_edge_len(double new_len){ m_max_edge_len = new_len;}

void set_texture_param(double surfTol, double angleTol, double max_edge_len=-1.){
	m_texture=true; m_tex_surfTol=surfTol;
	m_tex_angleTol=angleTol;m_tex_max_edge_len=max_edge_len;
}
bool texture()const{return m_texture;};
double tex_surfTol()const{return m_tex_surfTol;};
double tex_angleTol()const{return m_tex_angleTol;};
double tex_max_edge_len()const{return m_tex_max_edge_len;};

private:
	double m_crvTol;	///<バウンダリのトレランス
	double m_surfTol;	///<平面とみなすトレランス
	double m_max_ratio;	///<最大アスペクト比
	MGCL::fan_kind m_fk;
		///< =SINGLE_TRIANGLE,
		///<	1 triangle/FAN(default) and STRIP for as many as posible triangles.
		///<	STRIP triangles may cover multiple rectangles.
		///< =MULTIPLE_TRIANGLES,
		///<	as many triangles as possible/FAN and STRIP for as many as posible triangles.
		///<	STRIP triangles may cover multiple rectangles.
		///< =SINGLE_TRIANGLE_NO_STRIP,
		///<	SINGLE_TRIANGLE, but STRIP triangles cover only one tessellated rectagle.
		///< =MULTIPLE_TRIANGLES_NO_STRIP,
		///<	MULTIPLE_TRIANGLES, but STRIP triangles cover only one tessellated rectagle.
	size_t m_minimum_tri;	///<Specify minimum number of triangles.
	double m_max_edge_len;	///<when max_edge_len<=0, this means no limits on an edge length.

	bool m_texture;	///<true if texture setting is necessary, false if not.
		///<The following m_tex_xxxx member data are valid only when m_texture is true.
	double m_tex_surfTol;///<texture surface tolerance.
	double m_tex_angleTol;///<texture angle tolerance.
	double m_tex_max_edge_len;///<maximum width and height of the texture plane.
};

#endif
