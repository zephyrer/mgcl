/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _mgTLparameter_HH_
#define _mgTLparameter_HH_

#include <vector>
#include "Tl/TLisectsList.h"

class MGObject;
class MGFace;
class MGFSurface;
class MGSurface;
class mgTLparameter;
class mgTLPoints;
class mgTLInputParam;

///mgTLparameter is a proprietry class for Face tessellation.
///mgTLparameter holds necessary parameter data for face tessellation.
///In the constructor of mgTLparameter(const MGFace&, double, double),
///all the parameters are initialized.
///This parameter must be hold during the use of mgTLRects(m_rects member data).
///mgTLRects will reference mgTLparameter data.
class MGCLASS mgTLparameter{

public:

MGDECL friend std::ostream& operator<< (std::ostream& out, const mgTLparameter& para);

/////////Constructor///////////////
mgTLparameter():m_face(0), m_surface(0),m_points(0),m_edgPerim(0){;}
mgTLparameter(
	const MGFace& f,
	mgTLPoints* tlpoints,
	double crvTol,
	double surfTol,
	double max_ratio,
	double max_edge_len=-1.);
mgTLparameter(
	const MGSurface& f,
	mgTLPoints* tlpoints,
	double crvTol,
	double surfTol,
	double max_ratio,
	double max_edge_len=-1.
);
mgTLparameter(
	const MGFSurface& obj,
	mgTLPoints* tlpoints,
	double crvTol,
	double surfTol,
	double max_ratio,
	double max_edge_len=-1.,
	bool textue=false,
	double tex_surfTol=-1.,
	double tex_angleTol=-1.,
	double tex_max_edge_len=-1.
);
mgTLparameter(
	const MGFSurface& obj,	///<テセレーションするフェイス
							///<Must be MGFace or MGSurface.
	mgTLPoints* tlpoints,
	const mgTLInputParam& param///<parameter for the tessellation.
);
mgTLparameter(const mgTLparameter& param2);

////////////Destructor///////////////
~mgTLparameter();

size_t edgPerimLen()const;
std::vector<bool>* get_edgPerim(){return m_edgPerim;};
const MGFace& get_face()const{return *m_face;};
const MGSurface& get_surface()const{return *m_surface;};
double get_max_ratio()const{return m_max_ratio;};
double get_max_ratio_sqr()const{return m_max_ratio_sqr;};
double get_ratio_error()const{return m_error;};
double get_max_edge_len()const{return m_max_edge_len;};
mgTLisectsList& get_uList(){return m_uList;};
mgTLisectsList& get_vList(){return m_vList;};
double get_surface_error()const {return m_surface_error;};
double get_tess_crvError()const {return m_tess_crvError;};
double get_tess_srfError()const {return m_tess_srfError;};
void get_texture_parameter(
	double& tex_surfTol, double& tex_angleTol, double& tex_max_edge_len
)const;
double get_UError()const {return m_uError;};
double get_VError()const {return m_vError;};
bool is_face()const{ return m_face!=0;};
double isect_uerror()const{return m_puerror;};
double isect_verror()const{return m_pverror;};
const double* suf_param_range()const{return m_surf_prange;};
mgTLPoints* tlpoints()const{return m_points;};
bool texture()const{return m_texture;};

/////////Operator oveload/////////

///Assignment.
mgTLparameter& operator=(const mgTLparameter&);

private:

	const MGFace* m_face;	///<Original face to tessellate.
	const MGSurface* m_surface;///<Original surface to tessellate.
							///<If m_face!=0, m_surface=m_face->surface();
	mgTLPoints* m_points;	///<(u,v) points holder
	std::vector<bool>* m_edgPerim;
		///<Array of number of perimeter loops, m_edgPerim[j][i].
		///<If an outer loop exists, the array size is one.
		///<m_edgPerim[j][i] indicates if j-th loop's i-th edge is a perimeter
		///<edge or not, for 0<=j<edgPerimLen(). That is, if on a perimeter or not. 
		///<If true, it means the edge is a perimeter edge.
	mgTLisectsList m_uList;
	mgTLisectsList m_vList;
		///<Intersection points vector(mgTLisects) for u=min and u=max(m_uList),
		///<or v=min and v=max(m_vList) will be set initially.
		///<During the construction of mgTLRects, other necessary data will be
		///<added by the constructor, and will be referenced to process m_rects
		///<data of the mgTLRects.
	double m_puerror, m_pverror;///<Parameter error used for intersection computation.
	double m_surface_error;	///<error square of m_face->surface()->parameter_error();
		///<parameter of mgTLRect::normalize();
	double m_surf_prange[4];///<Surface parameter range, ={umin,umax,vmin,vmax}.
	double m_tess_crvError;	///<Tessellation curve tolerance.
	double m_tess_srfError;	///<Tessellation surface tolerance.
	double m_uError, m_vError;///<Parameter range error along u or v.
		///<Subdivision will never be performed that makes the result parameter
		///<range become less than above error. This is juded in mgTLdivideU().
	double m_error;			///<Ratio error.
		///<Subdivision to adjust ratio(ratio_subdivide) will never be performed
		///<when length of the rectangle along u or v is smaller than above m_error.
	double m_max_ratio;		///<Maximum ratio of rectangles' u and v sapan.
	double m_max_ratio_sqr;	///<Square of m_max_ratio.
	double m_max_edge_len;	///<Minimum edge length of the rectangles.

	bool m_texture;	///<true if texture setting is necessary, false if not.
		///<The following m_tex_xxxx member data are valid only when m_texture is true.
	double m_tex_surfTol;///<texture surface tolerance.
	double m_tex_angleTol;///<texture angle tolerance.
	double m_tex_max_edge_len;///<maximum width and height of the texture plane.
};

#endif
