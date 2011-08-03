/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// MGConstructionPlane.h : MGConstructionPlane クラスの宣言およびインターフェイスの定義をします。
#ifndef _MGConstructionPlane_HH_
#define _MGConstructionPlane_HH_

class MGBox;
class MGPosition;

#include "mg/Plane.h"
#include "mgGL/Color.h"

/** @addtogroup DisplayHandling
 *  @{
 */

///MGConstructionPlane defines a construction plane,
///which provides a local working 2D coordinate system and a local 3D coordinate system.
///MGConstructionPlane has the right hand coordinate system (U,V,N), where
///U=uspan*m_plane.uderiv(), V=uspan*m_plane.vderiv(), and N=m_nspan*m_plane.normal().
///m_plane.uderiv(), m_plane.vderiv(), and m_plane.normal() are set to length one
///on the construction.
///MGConstructionPlane has the cplane coordinate system such that the coordinate (x, y, z) is
///the normal world coordinate (R+x*U, R+y*V, R+z*N). The cplane coordinate conversion
///utilities are provided as convert_to_world() or convert_from_world().
class MGCLASS MGConstructionPlane{

public:

//////////Constructor/////////

MGConstructionPlane():m_disabled(true){;};
MGConstructionPlane(
	double origin[3],	///<Origin's coordinate value.
	double uaxis[3],	///<A vector value of the horizontal direction.
	double vaxis[3],	///<A vector value of the vertical direction.
	double uspan,		///<span length along u axis.
	double vspan,		///<span length along v axis.
	int uline_num,		///<number of lines along u axis.
	int vline_num,		///<number of lines along v axis.
	double nspan=1.		///<span length along normal axis.
);
MGConstructionPlane(
	const MGPlane& plane,///<construction plane.
	double uspan,		///<span length along u axis.
	double vspan,		///<span length along v axis.
	int uline_num,		///<number of lines along u axis.
	int vline_num,		///<number of lines along v axis.
	double nspan=1.		///<span length along normal axis.
);

///Construct a construction plane froma a box of 3D space.
MGConstructionPlane(
	const MGBox& box,
	int view_num=1	///<Standard view number:
	///<1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	///<0: non standard view.
);

//Copy constructor.
//MGConstructionPlane(const MGConstructionPlane& pobj2);

////////Destructor///////////
//~MGConstructionPlane();

///////////Operator overload////////////

//Assignment operator.
//MGConstructionPlane& operator=(const MGConstructionPlane& pobj);

////////////オペレーション/////////

///Bind the input point uv(this construction plane's parameter value)
///to the nearest grid point of this construction plane.
void bind_to_grid(
	const MGPosition& uv,
	MGPosition& uvout
)const;

///Convert cplane coordinates to the normal world coordinates.
///Function's return is the world coordinates.
MGPosition convert_to_world(const MGPosition& cplane_coord)const;

///Convert to cplane coordinates from the normal world coordinates.
///Function's return is the cplane coordinates.
MGPosition convert_from_world(const MGPosition& world_coord)const;

///Draw this plane using OpenGL.
void draw()const;

///Return if locate point on this plane should be bind to grid point or not.
///true if should be bound to grid point.
bool is_bind_to_grid()const{return m_bind_to_grid;};

///Obtain the position data of the parameter (u,v).
MGVector eval(const MGPosition& uv)const{return m_plane.eval(uv);};
MGVector eval(double u, double v)const{return m_plane.eval(u,v);};

///Get line and axis colors
void get_color(
	MGColor& lineColor,		///<Grid line color
	MGColor& uaxisColor,	///<u axis
	MGColor& vaxisColor		///<v axis
)const;

///locate a point on this plane, given straight line.
///the located point will be the intersection(or nearest bound) point
///of sl and the plane.
///Function's return value is the world coordinate located.
MGPosition locate(
	const MGStraight& sl,///<input the ray straight line of the cursor.
	MGPosition& uv		///<the plane's parameter value (u,v) will be output.
)const;

bool disabled()const{return m_disabled;};
bool enabled()const{return !m_disabled;};

///set bind_to_grid enable or disable.
void set_bind_to_grid_enable(){m_bind_to_grid=true;};
void set_bind_to_grid_disable(){m_bind_to_grid=false;};

///Compute grid data and the plane from box and view.
void set_grid_data(
	const MGBox& box,
	int view_num	///<Standard view number:
	///<1: 3D perspective view, 2:(x,y) 2D view, 3:(y,z) 2D view, 4:(z,x) 2D view
	///<0: non standard view.
);

///Compute grid data and the plane from the plane and the grid span data.
void set_grid_data(
	const MGPlane& plane,///<construction plane.
	double uspan,		///<span length along u axis.
	double vspan,		///<span length along v axis.
	int uline_num,		///<number of lines along u axis.
	int vline_num,		///<number of lines along v axis.
	double nspan=1		///<span length along normal axis.
);

///retrun the gl display list name of this construction plane.
size_t glname()const{return (size_t)&m_plane;};

bool valid()const{return m_plane.sdim()>0;};
const MGPlane& plane()const{return m_plane;};
MGPlane& plane(){return m_plane;};
double uspan()const{return m_uspan;};
double vspan()const{return m_vspan;};
int vnum()const{return m_vnum;};
int unum()const{return m_unum;};

///set the colors to default ones.
void set_default_color();

///set grid line color.
void set_line_color(const MGColor& color);

///set uaxis color.
void set_uaxis_color(const MGColor& color);

///set vaxis color.
void set_vaxis_color(const MGColor& color);

void set_disable(){m_disabled=true;};
void set_enable(){m_disabled=false;};
void set_span(double span){m_uspan=m_vspan=m_nspan=span;};
void set_uspan(double span){m_uspan=span;};
void set_vspan(double span){m_vspan=span;};
void set_num(int line_num){m_vnum=m_unum=line_num;};
void set_unum(int unum){m_unum=unum;};
void set_vnum(int vnum){m_vnum=vnum;};

private:
	bool m_disabled;///<true if this construction plane should not be displayed.
					///<false if this construction plane should be displayed.
	bool m_bind_to_grid;///<true if locate point should be bind to grid point of this plane.
					///<false if not.
	MGPlane m_plane;///<construction plane expression.
	double m_vspan;///<the span length of the direction m_plane.v_deriv().
	double m_uspan;///<the span length of the direction m_plane.u_deriv().
	double m_nspan;///<the span length of the direction m_plane.normal();
	int m_vnum;	///<number of lines along v axis.
	int m_unum;	///<number of lines along u axis
	MGColor m_lineColor,m_uaxisColor,m_vaxisColor;///<color of grid lines, u-axis and v-axis, each.

///Compute cplane parameter from 3D box.
MGDECL friend void MGcplane_parameter(
	const MGBox& box,
	double& span,	///<span length will be output.
	size_t& lnum,	///<number of lines along vertical and horizontal will be output.
	size_t& sdid,	///<maxmum area coordinate pair will be output.
					///<0:(x,y), 1:(y,z), 2:(z,x)
	MGPosition& mid	///<rounded mid point will be output.
);

///Debug Function.
MGDECL friend std::ostream& operator<< (std::ostream& out, const MGConstructionPlane& pln);

/// Serialization.
MGDECL friend MGOfstream& operator<< (MGOfstream& buf, const MGConstructionPlane& cpl);
MGDECL friend MGIfstream& operator>> (MGIfstream& buf, MGConstructionPlane& cpl);

};

/** @} */ // end of DisplayHandling group
#endif
