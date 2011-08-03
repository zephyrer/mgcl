/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLRect_HH_
#define _mgTLRect_HH_

#include "mg/Plane.h"
#include "Tl/TLRecisects.h"
#include "Tl/TLTriangles.h"

class mgTLPoints;
class mgTLRects;
class mgTLparameter;
class MGFace;
class mgTLData;

#if defined(MGCL_DLL)
#include "mg/ListProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGListProxy<mgTLRect*>;
#pragma warning( pop )
#else
#include <list>
#endif
#include <deque>

///mgTLRect is a proprietry class for Face tessellation.
///mgTLRect holds all the necessary information for the triangulation
///of a retangle face parameter space.
class MGCLASS mgTLRect{

public:

	///m_Recisects are intersections of the rectangle with a loop.
	///In m_Recisects, intersections are sorted in the anti-clock-wise along
	///the rectangle perimeter.
	///On return from mgTLRects constructor,
	///m_Recisects[2*j] is the point that is going into the rectangle and
	///m_Recisects[2*j+1] is going out. These two points always makes a pair.
	///(This is not the fact during the construction process of mgTLRects.
	///mgTLRect::normalize_trim_points in the constructor will take care
	///the sequence.)
	mgTLRecisects m_Recisects;	

#if defined(MGCL_DLL)
	typedef MGListProxy<mgTLRect*> container_type;
#else
	typedef std::list<mgTLRect*> container_type;
#endif

	///Iterator for m_neighbours.
	typedef container_type::iterator RecNItr;
	typedef container_type::const_iterator CRecNItr;

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLRect& rect);
MGDECL friend void change_neighbours(
	mgTLRect* rect2,
	size_t peri,
	std::list<mgTLRect*>& neibrs1,
	mgTLRect::RecNItr N1Itr);
//////////// constructor ///////////////
mgTLRect();

///Construct from the parameter box (u0, v0) - (u1, v1)
mgTLRect(double u0, double u1, double v0, double v1);

///Construct a higher parameter value mgTLRect by subdividing the input rect
/// along u(uDirection=ture) or v direction. This constructor builds
///the following data:
///m_urange, m_vrange, m_Pid, m_Pin, m_status.
///m_Pin[0], [3] (when uDirection=true) or m_Pin[0], [1] are not set.
mgTLRect(bool uDirection, mgTLRect& rect);

//////////// Member Function ///////////////

///Add vertices of the trimming line(boundary line) from isects()[sztRectIndex] to
/// [sztRectIndex+1], converting the boundary line to a polyline.
void add_boundary_vertices(
	double crvTol,
	size_t iss,
	mgTLTriangle& polygon,
	mgTLPoints& Tlpoints
)const;

///Add vertices of the neighboring rectangles' vertices at the perimeter peri
///from parameter s1 to s2. s1 and s2 are in the order of the rectangle. That is
///if peri==0 or 1, in the order of the u or v parameter value, and
///if peri==2 or 3, in the reverse order or the u or v parameter value.
void add_neighbor_vertices(
	size_t peri,
	double s1,
	double s2,
	mgTLTriangle& polygon
)const;

///add (u,v) to by add_point of rects.
size_t add_rect_point(
	size_t pnum,		///<perimeter number of this rect.
	double u,
	double v,
	mgTLRects& rects	///<mgTLRects into which all the rects are stored.
);

///get center(surface parameter) of this rect.
MGPosition center_uv()const;

///compare which is greater, t1 or t2 at the perimeter perim.
///latter position around the anti-clockwise round of the rectangle is greater.
///If t1 is greater or equal to t2 return true.
bool compare(size_t perim,double t1,double t2)const;

///Compute the closeset plane of this rect.
MGPlane compute_plane(const MGSurface& surf)const;

///Compute texture coordinates of the point world coord point xyz.
///This is valid only after texutured.
MGPosition compute_perimeter_tex_coord(
	const mgTLData& tldata,
	size_t i,///<perimeter number where uv lies on.
	const MGPosition uv	///<parameter value to get the texture, which
			///<must be on the perimeter i.
)const;

///Assurme this is not textured yet, and at least one of the
///neighboring perimeter is textured. Then compute_texture_by_triangles
///will perform texture computation by invoking mgTriangle::texture.
void compute_texture_by_triangles(
	mgTLData& tldata
);

///neighbor is already textured, and neighbor  and this is connected
///along neghbor's perimeter number nperim. Then
///compute_texture_from_neighbor will texture-compute this rect.
void compute_texture_from_neighbor(
	mgTLData& tldata,
	mgTLRect& neighbor,
	size_t nperim
);

///Assume four corner vertices are textured, texture other vertices.
void compute_texture_of_non_corner(
	mgTLData& tldata
);

///矩形交点からトリムポリゴンを作成し polygons に追加する
///トリムされている場合の処理を行う。頂点は mgTLPoints に追加する
void createPolygon(
	mgTLparameter& param,
	mgTLTriangles& polygons
);

///Get the end point parameter value of the perimeter peri.
///End point means in the order of the rectangle' anti-cloclwise order.
double end_point(size_t peri)const;

///Find the 1st neighbour iterator i whose minimum value at the perimeter peri
///is less or equal to t and the maximum value is greater than t.
///Function's return value is true if t is equal to minimum value of the perimeter.
bool find_neighbour(size_t peri, double t, RecNItr& i);

///Find the neighbor rect or face of rect_next following the line of fpoints.
///Function's return value is id of hhis to process the next face(mgTLData).
///When hhis moves to a different face, rect_next=0 and the proper next hhis id will
///be output.
int find_neighbor_rect(
	double errorSurf,///<surface parameter space error.
	int id_hhis,	///<id of hhis that defines the current priority line of uv is input.
					///<When id_hhis>=hhis.size(), no priority lines are input.
	const MGHHisect& hhis,
	int perim_into,	///<The perimeter number(of this) that came into this rect.
					///<perim_into=-1 when no previous rect.
	mgTLRect*& rect_next,///<next rect will be output.
						///<When rect_next=0 is output, no next in this face.
	int& perim_going_out///<rect_next's perimeter number.
);

///Get the point id of uv.
///get_point_id will() tests if uv is equal to a corner data of this rect, and
///if the same, the corner point id will be returned, else uv will be
///added to Tlpoints.
size_t get_point_id(
	const MGPosition& uv,	///<parameter value to get the texture, which
			///<may not be on the perimeter i.
	mgTLPoints& Tlpoints
)const;

///Compute texture coordinates of the point world coord point xyz.
///This is valid only after texutured.
MGPosition get_tex_coord(
	const mgTLData& tldata,
	const MGPosition& uv	///<parameter value to get the texture, which
			///<may not be on the perimeter i.
)const;

///test if the parameter uv is on a perimeter of this rect.
///returned is:
///=-1: when not on any perimeter,
///>=0: when on a perimeter, perimeter id will be returned.
int is_on_perimeter(const MGPosition& uv)const;

///test if the parameter uv is on a corner of this rect.
///returned is:
///=-1: when not on any corner,
///>=0: when on a corner, the corner id will be returned.
int get_corner_id(const MGPosition& uv)const;

///Test if finish triangulation, return true.
bool is_triangled()const{return m_triangle!=-1;}

///Get the starting triangle id.
int triangle_id()const{return m_triangle;}

///rectがストリップになる条件を持っているかどうか調べる
///条件 ・隣り合うrectが1以下 ・トリムされていない ・triangulationされていない
bool hasStripCondition()const;

///Test if this rect has trimming intersection or not(for the 4 perimeters).
bool has_trim(){return m_Recisects.size()>1;}

///Test if this rect has trimming intersection or not on the perimeter peri.
bool has_trim(size_t peri) const;

///隣のrectの頂点を求める。頂点はパラメータの小さい順に並んでいる。
void getNeighborPoint(size_t peri, std::vector<size_t>& vecPoint) const;

///Test if this rectangle includes the parameter value uv.
///Returns true if includes uv, false if not.
bool includes(const MGPosition& uv)const;

///This is another version of includes, and takes care of tolerance.
bool includes2(const MGPosition& uv)const;

///Make the uv included in this rectangle's parameter space.
void make_included(MGPosition& uv)const;

///ペリメタ番号とパラメータ値を与えて次のスタート交点を求める
///求まったときtrueが返却されindexに交点インデックスが入る
bool NextStartIsect(
	size_t peri,///<perimeter number to get.
	double p,	///<last isect's parameter value.
	size_t periS,///<polygon's starting perimeter number and 
	double sS,	///<the parameter value.
	size_t& index
)const;

///Compute trimming points(intersections with trimming loops with the rectangle)
///at the newly generated adjacent perimeters when old rectangle is subdivided
///into half at the parameter t(u or v according to the coordinate kind kcod).
///Trimming points will be recomputed and stored in this and rect2's m_Recisects.
void get_trim_point(
	mgTLparameter& param,	///<Tessellation parameter
	double t,				///<Subdividing parameter value u or v.
	size_t kcod,			///<t's coordinate kind. 0:u, 1:v.
	double s0, double s1,	///<parameter range of the other coordinate than t.
	mgTLRect& rect2			///<2nd rectangle after subdivided.
);

///Return i-th perimeter's neighbours.
const std::list<mgTLRect*>& neighbours(size_t i)const{
#if defined(MGCL_DLL)
	return m_neighbours[i].std_container();
#else
	return m_neighbours[i];
#endif
}

size_t opposite_perim(size_t i){ return (i+2)%4;};

size_t Pid(size_t i)const{return m_Pid[i];};
const size_t* Pid()const{return m_Pid;};

bool Pin(size_t i)const{return m_Pin[i];};
const bool* Pin()const{return m_Pin;};

void print_neighbours(std::ostream& out) const;
void print_triangles(std::ostream& out, mgTLData& tldata) const;

///Compute ratio square of du/dv of the rectangle.
double ratio_sqr(const mgTLparameter& param)const;

///Obtain previous or aft ReciItr of i in the ReciItr loop of this
///rectangle.
mgTLRecisects::ReciItr reci_aft(mgTLRecisects::ReciItr i);
mgTLRecisects::ReciItr reci_pre(mgTLRecisects::ReciItr i);

///Set all the texture coordinates data of this rect's vertex.
void set_initial_tex_from_1point(
	const MGPosition& xyz1,	///<parameter value of the surfaceof texture st1.
	const MGPosition& st1,	///<texture coordinates of uv1
	const MGUnit_vector& xaxis,	///<texture x-axis.
	mgTLData& tldata
);

///Set the new m_status of mgTLRect after subdivided.
///m_box, m_Rectisects, and m_Pin[0] must be updated before use.
///These are referenced in set_status(). old m_status is referenced also.
void set_status(
	const MGFace& f	///<Original face.
	);

///Assuming this is already textured, 
///set all the neighboring rect's texture coordinates data
///invoking compute_texture_by_triangles().
///Function's return value is true if some are textured.
bool set_neighbor_tex_coord_data(
	mgTLData& tldata,
	std::deque<mgTLRect*>& rects///<modified rects will be prepended.
);

///set triangled flag.
void set_triangled(int tri);

///Get the start point parameter value of the perimeter peri.
///Start point means in the order of the rectangle' anti-cloclwise order.
double start_point(size_t peri)const;

///Return the status of the rectangle.
int status() const{return m_status;};

///Return the m_Recisects reference.
const mgTLRecisects& isects() const{return m_Recisects;}
mgTLRecisects& isects(){return m_Recisects;}

///Normalize trimming points.
///Re-order the recisect seq 
void normalize_trim_points();

///Divide into two parts at the middle of u(uDirection=true) or
///v(uDirection=false) parameter.
///The divided two are put in this and function's return value.
///The new this will be located first in the sequence of
///the uv rectangle list mgTLRects.
mgTLRect* subdivide(
	mgTLRects& rects,	///<mgTLRects into which all the rects are stored.
	bool uDirection	///<direction of sibdivision. =true:udirection.
	);

///Test if texture coord data is already set.
bool is_textured()const{return m_textured>0;}

///Return this rectangle's corner point i's (u,v).
MGPosition uv(size_t i)const;

///Return corner point (u,v) data.
double umin()const{return m_u0;};
double umax()const{return m_u1;};
double uspan()const{return m_u1-m_u0;}
double vmin()const{return m_v0;};
double vmax()const{return m_v1;};
double vspan()const{return m_v1-m_v0;}

void set_tex_within_tol(bool within_tol=true){m_tex_within_tol=within_tol;};
bool tex_within_tol()const{return m_tex_within_tol;};
void set_texture_level(int level){m_textured=level;};
int texture_level(){return m_textured;};

private:
	int m_status;	///<Status of the rectangle: in, out, on, over, on&over.
		///<	1:MGRECT_IN,	//Whole rectangle is inside the Face.
		///<	2:MGRECT_OUT,	//Whole rectangle is outside the Face.
		///<	3:MGRECT_ON,	//Some trimming curve is crossing the rectangle.
		///<					//So part is inside, and part is outside the face.
		///<	4:MGRECT_OVER,	//In the case that whole inner boundary loop is included
							///<in the rectangle, perimeter loops exists, or at least
							///<one edge of the outer loop is not on the surface perimeter.
		///<	5:MGRECT_ONANDOVER
			///<In the rectangle whole inner boundary loop is included, and
			///<Some trimming curve is crossing the rectangle.
	
	///Number of subdivision.
	int m_divnum;///<Initially 0, and whenever the subdivision occured,
				///<will be incremented.

	double m_u0,m_u1,m_v0, m_v1;///<Parameter range of the rectangle, that is,
			///<m_u0-m_u1:u range, m_v0-m_v1:v range.

///The following variables' array id is a vertex or perimter number of the rectangle.
	bool m_Pin[4];		///<m_Pin[i] indicates if the position i is in the range of
						///<the face or not. The i is a vertex number of the rectangle.
	size_t m_Pid[4];	///<m_Pid[i] is the index of mgTLRects' m_points.
						///<(*m_points)[m_Pid[i]] is the i-th vertex postion data.
	container_type m_neighbours[4];
		///<m_neighbour[i] is a list of neighbour rects of this rect's perimeter i.
		///<In m_neighbour[i], neighbor rects will be sorted in the order of
		///<the perimeter's parameter value. That is, if i=0 or 2, in the order of
		///<the parameter value of u. If i=1 or 3, in the order of the parameter v.

	int m_triangle;///<After triangulated, id of tldata's m_triangles will be set.
		///<if not, m_triangle=-1.

	int m_textured;	///<=0:if not textured yet at all.
					///<=n:level n textured. Level 1 means texture processed from fix line
					///<		or a point.
					///<	 Level n means texture processed from neighbor level n-1 rect.
					///<   Here n>=2.
	double m_texcoord_ratio;///<ratio of world length and texture coordinate length.
		///<initially 1. is set.
	bool m_tex_within_tol;///<true if this rect is already within tolerance of a plane for
		///<texture mapping.

///Add rect2 after rect1 in the neibour sequence of this perimeter peri.
void add_neighbour(
	size_t peri,		///<Perimeter num of this.
	mgTLRect* rect1,
	mgTLRect* rect2);

///Test if all corner points are in the face or not.
///Return true if all corner points are in.
bool all_corners_in()const;

///Change neighbours of this rect's neibrs1 from N1Itr to neibrs1.end() to neibrs2's
///neighbours. Here, neibrs2 is rect2->m_neighbours[peri], and
///neibrs1 is this->m_neighbours[peri].
void change_neighbours(
	mgTLRect* rect2,
	size_t peri,
	mgTLRect::RecNItr N1Itr
);

///Test if this rect has perimeter boundary, check the perimeter boundary
///is within the curve tolerance error.
bool checkPerimeterBoundary(
	mgTLPoints* points,
	double crvTol,
	const MGSurface &srf,
	bool &bSubDivU);

///矩形からポリゴンを作成し polygon を作成し、polygons に追加する
///トリムされていない場合の処理を行う
void createNonTrimPolygon(mgTLTriangles& polygons);

///矩形交点からトリムポリゴンを作成し polygons に追加する
///トリムされている場合の処理を行う。頂点は mgTLPoints に追加する
void createTrimPolygon(
	mgTLparameter& param,
	mgTLTriangles& polygons
);

///divide neighbours at the parameter value t of the perimeter num peri,
///and put the latter part to rect2.
void divide_neighbours(size_t peri, double t, mgTLRect& rect2);

///Check if this rect's perimeter(that is rect's v=min,max, when along_u is true
///or u=min,max parameter line, when aling_u is false) is within error .
bool perim_within_tol(
	bool along_u,
	const MGSurface& srf,
	double error
)const;

bool valid()const;

friend class mgTLRects;
friend class mgTLData;

};

#endif
