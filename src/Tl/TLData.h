/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLData_HH_
#define _mgTLData_HH_

class MGObject;
class MGFSurface;
class MGSurface;
class MGFace;
class mgTLInputParam;
class mgTLRect;
class mgTLTexPoints;
class mgTLTriangle;
class mgTLTriangles;
class mgTLDataVector;
class mgTLTexPlane;
class MGTextureImage;
class MGTextureImage2;
class MGTextureImages;

#include <deque>
#include "Tl/TLparameter.h"
#include "Tl/TLRects.h"
#include "Tl/TLPoints.h"
#include "Tl/TLTexPoints.h"
#include "Tl/TLTriangles.h"

/** @addtogroup UseTessellation
 *  @{
 */

/// A class that contains one Face or Surface's tessellated data.
///mgTLData includes the following two newed pointers:
///	mgTLparameter* m_tlparam;	//newed tl parameter.
///	mgTLRects* m_rects;			//newed tl rects.
///	mgTLTriangles* m_triangles;	//newed mgTLTriangles pointer.
///	mgTLPoints* m_tlpoints;		//newed mgTLPoints pointer.
///These pointers are treated as std::auto_ptr. When copied, or assigned, the ownership
///are transfered to the new mgTLData object and the original object's pointers are set
///null.
class MGCLASS mgTLData{

public:

enum TEXTURE_STATUS{
	NOT_TEXTURED=0,
	PARTIALLY_TEXTURED,
	TOTALLY_TEXTURED
};

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLData& data);

/////////////Constructor////////////
mgTLData(const mgTLData& tlData);///copy constructor.

///Construct by tessellating the obj.
///To rendering by shading model, use mgGDL::shade().
mgTLData(
	const MGFSurface& obj,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	const mgTLInputParam& param///<parameter for the tessellation.
);

///Construct by tessellating the obj.
///To rendering by shading model, use mgGDL::shade().
mgTLData(
	const MGFSurface& obj,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	double crvTol,			///<バウンダリのトレランス
	double surfTol,			///<平面とみなすトレランス
	double max_ratio=2.,	///<最大アスペクト比
	MGCL::fan_kind fk=MGCL::MULTIPLE_TRIANGLES,
		///<fk=SINGLE_TRIANGLE:   1 triangle/FAN
		///<fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri=8,	///<Specify minimum number of triangles.
	double max_edge_len=-1.	///<when max_edge_len<=0, this means no limits on an edge length.
);

///////////// Destructor./////////
~mgTLData();

///////////// Operator overload./////////

/// Assignment
mgTLData& operator=(const mgTLData& tld);

#if defined(MGCL_DLL)
bool operator<(const mgTLData&) const{ return true;}
bool operator==(const mgTLData&) const{ return true;}
#endif /// MGCL_DLL

////////// member function /////////////

///compute all the texture coordinates of this m_tlpoints.
///compute_texture_coordinates can be invoked only after triangulate().
void compute_texture_coordinates(
	 const MGPosition& uv,	///<surface parameter of the point st.
	 const MGPosition& st,	///<texture coordinate of uv.
	 const MGUnit_vector& texXaxis,///<world coordinate x axis vector.
 	 double distortion_ratio=-1.///<When >1. distorted points(points whose ratio are more
							///<than distortion_ratio or less than 1/distortion_ratio will be
							///<stored in m_tex_distorted.
);

///compute the texture coordinates of the point uv from the data of 2 points,
///(uv1, st1) and (uv2, st2). Here uvx means the surface parameter and stx
///does texture coordinates known.
void compute_uv_texture_coordinate(
	const MGVector& normal,	///<Normal closest plane of the objective mgTLRect.
	 const MGPosition& uv,	///<surface parameter of the point st to compute.
	 const MGPosition& uv1,	///<surface parameter of the point st1.
	 const MGPosition& st1,	///<texture coordinate of uv1.
	 const MGPosition& uv2,	///<surface parameter of the point st2.
	 const MGPosition& st2,	///<texture coordinate of uv2.
	 MGPosition& st	///<texture coordinate of uv will be output.
)const;

///compute all the texture coordinates of this m_tlpoints.
///compute_texture_coordinates can be invoked only after triangulate().
///Function's return value is true if textured.
bool compute_texture_coordinates_partially(
	int& id_hhis,	///<id of hhis that define the current priority line of uv is input.
					///<When id_hhis>=hhis.size(), no priority lines are input, and
					///<all the rects will be textured from the point uv.
					///<id of hhis to process next face(mgTLData) will be output.
	const MGHHisect& hhis,
	const MGPosition& uv,
	const MGPosition& st,	///<texture coordinate of uv.
	const MGUnit_vector& texXaxis,///<world coordinate x axis vector.
 	double distortion_ratio=-1.///<When >1. distorted points(points whose ratio are more
							///<than distortion_ratio or less than 1/distortion_ratio will be
							///<stored in m_tex_distorted.
);

///compute all the non corner texture coordinates of this m_tlpoints.
///compute_non_corner_texture_coordinates can be invoked only after
///all the rects are textured.
void compute_non_corner_texture_coordinates();

///compute texture coordinates, given parameter edge connected to already
///computed face.
void compute_edge_texture_coordinates(
	const MGEdge& this_edge,	///<edge of this TLData connected to already computed face.
	const MGEdge& partnerE,		///<a partner edge of this_edge. The partner edge is already
		///<computed face's edge.
	mgTLData& partner_tld		///<partnerE's TLData which is already texture-computed.
);

///true when texture distortion record is required.
bool distortion_record_required()const{return m_tex_distorted!=0;};

std::vector<MGPosition>* distortion_record_vector(){return m_tex_distorted;};

///Draw textures within the m_texture_planes[plane_id].
void drawt_texture_image_part(
	const MGTextureImage2& teximage,	///<texture image data.
	const double back_color[4],	///<back ground color data.
	int plane_id	///<indicates if only parts of texture data be drawn.
		///<let i1=m_texture_planes[plane_id]->index(), and i2=m_texture_planes[plane_id+1]->index().
		///<Then rects of m_rects[i1] to m_rects[i2-1] are to be drawn.
)const;

///Draw textures within the m_texture_planes[plane_id].
void drawt_texture_image_part(
	const MGTextureImages& teximage,	///<texture image data.
	const double back_color[4],	///<back ground color data.
	int plane_id	///<indicates if only parts of texture data be drawn.
		///<let i1=m_texture_planes[plane_id]->index(), and i2=m_texture_planes[plane_id+1]->index().
		///<Then rects of m_rects[i1] to m_rects[i2-1] are to be drawn.
)const;

///find the rect that includes the parameter uv.
mgTLRect* find_rect(const MGPosition& uv)const{return m_rects->find_rect(uv);};

///Assumed this is already textured, compute texture coordinate at the point uv,
///and the world coordinate direction of the s-axtis.
///Function's return value is true if data is obtained successfully,
///false if the rect at uv is not textured and data was not obtained.
bool get_texXaxis(
	const MGPosition& uv,	///<parameter value of the face to compute at.
	MGPosition& st,			///<texture coordinate will be output.
	MGUnit_vector& texXaxis2///<the s-axis direction at uv in world coordinate.
)const;

///Search neighboring faces that is alreday texture computed,
///then after contact rectagles texture computation, compute the whole recttangle
///in this tlData.
///This is assumed not to be texture-computed yet.
///Function's return value is true if neighboring texture-computed face is
///found and texture-computed.
///false if not texture-computed.
bool find_textured_and_compute(
	mgTLDataVector& tldvec,
 	double distortion_ratio///<When >1. distorted points(points whose ratio are more
							///<than distortion_ratio or less than 1/distortion_ratio will be
							///<stored in m_tex_distorted.
);

///Compute minimum and maximum curvature of MGCL::SURFACE_CURVATURE_KIND kind.
void get_min_max_curvatures(
	MGCL::SURFACE_CURVATURE_KIND kind,
	double& minimum,
	double& maximum
)const;

///Get texture coordinates, given m_tex_coords id i.
const MGPosition& get_texcoord(size_t i)const;
MGPosition& get_texcoord(size_t i);

///Get surface parameters, given m_tlpoints id i.
const MGPosition& get_uv(size_t i)const ;
MGPosition& get_uv(size_t i);

///Get texture coordinates and world coordinates, given m_tlpoints id i.
void get_texcoord_world(
	size_t i,
	MGPosition& st,	///<(s,t) will be output, which is biased by base.
	MGPosition& world,///<world (x,y,x) value of the point (s,t) will be output.
	const MGVector& base
)const;

///double image_width()const{return m_tex_coords->image_width();};
///double image_height()const{return m_tex_coords->image_height();};

///Obtain the number of triangles tessellated.
size_t number_of_triangles()const;

///Obtain the number of points of the i-th triangle.
size_t number_of_points(size_t i)const;

mgTLparameter& parameter(){return m_rects->parameter();};

///Record the distortion point.
///uvr is (uvr[0], uvr[1]) is the surface parameter,
///and uvr[2] is the distortion ratio.
void record_tex_distortion(const MGPosition& uvr);

///Obtain i-th rect reference.
const mgTLRect& rect(int i)const{return m_rects->rect(i);};

void set_texcoord(size_t i, double s, double t);
void set_texcoord(size_t i, const MGPosition &texc);
void set_texcoord_as_not_modify(size_t i, const MGPosition &texc);

///Set texture computation status.
void set_textured(TEXTURE_STATUS texture_status){m_textured=texture_status;};

///Obtain the MGFace pointer if this is the data of MGFace,
///else null will be returned.
const MGFace* face_pointer()const;

///Obtain the surface& of this rect.
const MGSurface& surface()const;

void triangulate(
	MGCL::fan_kind fkind		///<実行パラメータデフォルトは1Fan,1Triangle
);

///Compute uv's texture coordinates.
MGPosition tex_coord(const MGPosition& uv);

bool tex_distorted(double ratio){
	return ratio>m_tex_distortion_ratio;
}

const std::vector<mgTLTexPlane*>& tex_planes()const{
	return m_rects->texture_planes();
}
const mgTLTexPlane* tex_plane(int i)const{
	const std::vector<mgTLTexPlane*>& texpls=tex_planes();
	return texpls[i];
}

///Compute maximum  width and height out of all the tex_planes.
void tex_plane_maximum(double& width, double& height)const;

///Test if this tlData's texture coordinates are computed
///invoking compute_texture_coordinates() .
///True will be returned if invoked.
TEXTURE_STATUS textured()const{return m_textured;};

const mgTLparameter& tlparam()const{return *m_tlparam;};
mgTLparameter& tlparam(){return *m_tlparam;};

const mgTLPoints& tlpoints()const{return *m_tlpoints;};
mgTLPoints& tlpoints(){return *m_tlpoints;};

const mgTLTexPoints& tlTexPoints()const{return *m_tex_coords;};
mgTLTexPoints& tlTexPoints(){return *m_tex_coords;};

const mgTLRects& tlrects()const{return *m_rects;};
mgTLRects& tlrects(){return *m_rects;};

///Obtain i-th triangle in this mgTLData.
const mgTLTriangle& triangle(size_t i)const;
mgTLTriangle& triangle(size_t i);

///Obtain the i-th triangle's type.
///=mgTESTRIANG_FAN:fan, =mgTESTRIANG_STRIP:strip.
mgTESTRIANG triangle_type(size_t i)const;

///Obtain mgTLTriangles pointer.
const mgTLTriangles& triangles()const{return *m_triangles;};
mgTLTriangles& triangles(){return *m_triangles;};

///Obtain the i-th traiangle's j-th point's surface parameter value.
const MGPosition& uv(size_t i, size_t j)const;
MGPosition& uv(size_t i, size_t j);

private:
	mgTLparameter* m_tlparam;	///<newed tl parameter.
	mgTLRects* m_rects;			///<newed tl rects.
	mgTLTriangles* m_triangles;	///<newed mgTLTriangles pointer. vector of Ids of m_tlpoints.
	mgTLPoints* m_tlpoints;		///<newed mgTLPoints pointer. Surface parameter(u,v).

	mgTLTexPoints* m_tex_coords;///<newed mgTLTexPoints pointer. Texture coordinates(s,t).
		///<m_tex_coords[i] is the texture coordinates of m_tlpoints[i], if m_tex_coords!=null.
	TEXTURE_STATUS m_textured;

	std::vector<MGPosition>* m_tex_distorted;///<Texture coordinate distorted points in this mgTLData.
		///<Let Pi=(*m_tex_distorted)[i], then (Pi[0],Pi[1]) is (u,v) of the surface parameter,
		///<Pi[2] is the distortion ratio at (u,v).
	double m_tex_distortion_ratio;///<distortion ratio. must be >1.

///Find the rects of level number (level-1)(which are supposed to be textured
///already) and texture their not-textured neighbor rects.
///Function's return value is
///false if no rects of (level-1) which were not textured were found.
///true if some rects of (level-1) which had non_textured neighbor rects
///were found and textured.
bool compute_level_texture(
	int level	///<Level number to texture.
);

///compute all the texture coordinates of this m_tlpoints.
///compute_texture_coordinates can be invoked only after triangulate().
void compute_texture_coordinates_from_neighbor();

///多角形からファンを作成して、m_Trianglesに追加する
void createMultiFanSet(
	mgTLTriangles& polygons
);

///多角形からSingleファンを作成して、m_Trianglesに追加する
void createSingleFanSet(
	mgTLTriangles& polygons
);

///与えられたrectから開始してv=maxの辺よりストリップを作成する
///ストリップが求まったらm_trianglesに追加される
///Function's return value is true if pRect is triangulated,
///false if not.
bool createStrip(
	mgTLRect* pRect
);

/// OpenGL display for shading.
void draw_texture(
	const MGTextureImage& teximage
)const;

void draw_texture(
	const MGTextureImages& timages
)const;

void draw_texture(
	const MGTextureImages& timages,
	mgTLTriangles::const_triIterator first,
	mgTLTriangles::const_triIterator last
)const;

///Draw textures fo tri. tir's texture coordinates are stored in tlTexPoints() as their
///woordl coordinate length. draw_texture_image_tri() draws texture of tri converting
///the texture coordinates to [0,1] span using the box(of the image).
///box provides the texture's (left, bottom) coordinates and the length of
///the image's [0,1] span.
void draw_texture_image_tri(
	const mgTLTriangle& tri,
	const MGBox& box///<The box (left,bottom, right, top) of the image
)const;

void draw_texture_image_tri_with_trim(
	const MGBox& image_box,///<The image box to draw.
	const mgTLTriangle& triangle,
	const MGBox& trim_box///<The trimming box of the triangle.
)const;

///Initialize m_tex_coords and all the rects' m_center and m_normal.
void initialize_tex_coords(
	double distortion_ratio
);

void set_up_texture_data(
	bool from_not_modify,///<true if only tlTexPoints data set_as_not_modify is to use.
	mgTLRect*& rect	///<1st modified rect will be output.
);

friend class mgTLDataVector;

};

/** @} */ // end of UseTessellation group
#endif
