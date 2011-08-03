/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLDatas_HH_
#define _mgTLDatas_HH_

#include <iosfwd>
#include <vector>
#include "Tl/TLData.h"

class MGShell;
class MGGroup;
class mgTLInputParam;
class MGFPoint;
class MGFSurface;

/** @addtogroup UseTessellation
 *  @{
 */

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<mgTLData>;
#pragma warning( pop )
#endif

typedef std::vector<MGFPoint> MGFPoints;

///mgTLDataVector is a vector of mgTLData which is a tessellated data of one face or surface.
///If two (or more than two) faces share a binder edge, the faces are guaranteed to
///be tessallated from the same one straight line.
class MGCLASS mgTLDataVector{

public:

	typedef std::vector<mgTLData>::iterator dataIterator;
	typedef std::vector<mgTLData>::const_iterator const_dataIterator;
	typedef std::vector<mgTLData>::iterator iterator;
	typedef std::vector<mgTLData>::const_iterator const_iterator;
	std::vector<mgTLData> m_datas;	///mgTLData sequence.

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLDataVector& datas);

//////////// constructor ///////////////

///void constructor
mgTLDataVector(){};

///construct an mgTLDataVector from an object.
///To rendering by shading model, use mgGDL::shade() for the member mgTLData.
mgTLDataVector(
	const MGObject& obj,	///<obj must be a MGSurface, MGFace, or MGShell
	const mgTLInputParam& param///<parameter for the tessellation.
);
mgTLDataVector(
	const MGObject& obj,	///<obj must be a MGSurface, MGFace, or MGShell
	double crvTol,			///<バウンダリのトレランス
	double surfTol,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

///construct an mgTLDataVector from an object.
///To rendering by shading model, use mgGDL::shade() for the member mgTLData.
mgTLDataVector(
	const MGGroup& group,	///<From group, a MGSurface, MGFace, or MGShell will be extracted.
	const mgTLInputParam& param///<parameter for the tessellation.
);
mgTLDataVector(
	const MGGroup& group,	///<From group, a MGSurface, MGFace, or MGShell will be extracted.
	double crvTol,			///<バウンダリのトレランス
	double surfTol,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

///mgTLDataVector(const mgTLDataVector& datas2);	///Copy constructor.

////////// destructor ///////////////
///~mgTLDataVector();

//////////// operator overload ////////////
///mgTLDataVector& operator= (const mgTLDataVector& datas2);

///Return the i-th mgTLData in this mgTLDataVector.
const mgTLData& operator[](size_t i)const{return m_datas[i];};
mgTLData& operator[](size_t i){return m_datas[i];};

////////// member function /////////////

const_dataIterator begin()const{return m_datas.begin();}
dataIterator begin(){return m_datas.begin();}
const_dataIterator end()const{return m_datas.end();}
dataIterator end(){return m_datas.end();}

void clear(){m_datas.clear();};

///compute all the texture coordinates of this m_tlpoints.
///compute_texture_coordinates will find f in this vector,
///then compute texture coordinates of f.
///If f is a constituent face of a shell, all the textures of the faces of the shell
///will be  computed.
void compute_texture_coordinates(
	const MGFSurface& f,	///<face or surface of uv.
	const MGPosition& uv,	///<surface parameter of the point st.
	const MGPosition& st,	///<texture coordinate of uv.
	 const MGUnit_vector& texXaxis,///<world coordinate x axis vector.
 	 double distortion_ratio=-1.///<When >1. distorted points(points whose ratio are more
							///<than distortion_ratio or less than 1/distortion_ratio will be
							///<stored in m_tex_distorted.
);

///compute all the texture coordinates of this m_tlpoints.
///compute_texture_coordinates will find f in this vector,
///then compute texture coordinates of f.
///If f is a constituent face of a shell, all the textures of the faces of the shell
///will be  computed.
void compute_texture_coordinates(
	const std::vector<MGFPoints>& fpointsVec,
				///<vector of vector of MGFPoint for high priority line.
				///<fpoints.size() must be >=1, and fpoints[0].size() must be >=1.
	const MGFSurface& f,	///<face or surface of uv.
	const MGPosition& uv,	///<surface parameter of the point st.
	const MGPosition& st,	///<texture coordinate of uv.
	const MGUnit_vector& texXaxis,///<world coordinate x axis vector.
 	double distortion_ratio=-1.///<When >1. distorted points(points whose ratio are more
							///<than distortion_ratio or less than 1/distortion_ratio will be
							///<stored in m_tex_distorted.
);

void draw_texture(
	const MGTextureImage2& teximage,
	const double* back_color
)const;

void draw_texture(
	const MGTextureImages& teximages,
	const double* back_color
)const;

bool empty()const{return m_datas.empty();};

///Find mgTLData that is for the input face in this vector.
///If found the pointer will be returned,
///else null will be returned.
mgTLData* find_tldata(const MGFSurface* face);

///Compute minimum and maximum curvature of MGCL::SURFACE_CURVATURE_KIND kind.
void get_min_max_curvatures(
	MGCL::SURFACE_CURVATURE_KIND kind,
	double& minimum,
	double& maximum
)const;

void push_back(const mgTLData& tldata){m_datas.push_back(tldata);};

///construct mgTLData by the input parameters and push backt the data.
void push_back(
	const MGObject& obj,	///<obj must be a MGSurface, MGFace, or MGShell
	const mgTLInputParam& param///<parameter for the tessellation.
);
void push_back(
	const MGObject& obj,	///<obj must be a MGSurface, MGFace, or MGShell
	double crvTol,			///<バウンダリのトレランス
	double surfTol,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

///construct mgTLData by the input parameters and push backt the data.
void push_back(
	const MGGroup& group,	///<From group, a MGSurface, MGFace, or MGShell will be extracted.
	const mgTLInputParam& param///<parameter for the tessellation.
);
void push_back(
	const MGGroup& group,	///<From group, a MGSurface, MGFace, or MGShell will be extracted.
	double crvTol,			///<バウンダリのトレランス
	double surfTol,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

size_t size()const{return m_datas.size();};

///Obtain i-th mgTLData in this mgTLDataVector.
const mgTLData& tldata(size_t i)const{return m_datas[i];}
mgTLData& tldata(size_t i){return m_datas[i];}

///Compute maximum  width and height out of all the tex_planes.
void tex_plane_maximum(double& width, double& height)const;

private:

///This DataVector is a vector of members of a shell, and at least one of
///mgTLData(that is, at least one of the member face) is textured.
///Then compute_shell_texture_coordinates textures all of non-textured mgTLData.
void compute_shell_texture_coordinates(
	double distortion_ratio
);

friend class mgTLData;

};

/** @} */ // end of UseTessellation group
#endif
