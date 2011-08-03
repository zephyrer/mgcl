/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#if defined(MGCL_NO_MFC)
#include <vector>
#include "mg/Tolerance.h"
#include "mg/DrawFunc.h"
#include "mg/Object.h"
#include "mg/Box.h"
#include "mg/LBRep.h"
#include "mg/SPointSeq.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Plane.h"
#include "mg/Group.h"
#include "mg/MGStl.h"
#include "topo/Face.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Shell.h"
#include "mgGL/Color.h"

class GdiplusStartupInput;

namespace mgGDL{
struct vertex3dv{void operator()(const MGPosition& pos) const{;};};
void setup_easy_shade(){}
void setup_wire_drawing(){}
void MGDrawArrow(MGPosition pos[4]){}
void MGDrawArrow(const MGCurve& curve, int ndiv){}
void MGDrawArrow(const MGFSurface& face, int udiv, int vdiv){}
void MGDrawBox(const MGBox& box){}
void MGDrawPointSeq(
	const MGBPointSeq& bp,
	bool draw_points,		//True if points be drawn.
	bool dotted				//true if dotted lines be drawn.
){}
void MGDrawPointSeq(
	const MGSPointSeq& sp,
	bool draw_points,		//True if points be drawn.
	bool dotted				//true if dotted lines be drawn.
){}
void MGDrawCurvaGraph(const MGCurve& crv, double scale, int density, bool use_radius){}
void MGDrawPoint(double x,double y,double z){}
void MGDrawPoint(const MGPosition& pos){}
void MGGLBeginLINE_STRIP(){;}
void MGGLVertex(double x, double y, double z){;}
void MGGLEnd(){;}
void param_rectangle(
	const MGBox& param_box	//Parameter range.
){}
void mgGDL_WTess(
	const mgTLData& tld	//tessellation data.
){}
void mgGDL_PTess(
	const mgTLData& tld	//tessellation data.
){}
void MGDrawPolyline(const MGBPointSeq& line, bool closed){}
void MGDrawPolyline(const std::vector<MGPosition>& line, bool closed){}
void MGDrawStraight(const MGPosition& end, const MGPosition& start){}
void draw_surface_curvature_data(
	const mgTLData& tld,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper, //minimum and maximum value of the curvatures of the kind.
		//Whne lower>=upper, lower is set as the minimum value and upper is set
		//as the maximum value out of all the curvatures.
	double* real_lower,	//When not null, minimum curvature of the kind will be output.
	double* real_upper	//When not null, maximum curvature of the kind will be output.
){}
void draw_surface_curvature(
	const mgTLDataVector& tldvec,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper, //minimum and maximum value of the curvatures of the kind.
		//Whne lower>=upper, lower is set as the minimum value and upper is set
		//as the maximum value out of all the curvatures.
	double* real_lower,	//When not null, minimum curvature of the kind will be output.
	double* real_upper	//When not null, maximum curvature of the kind will be output.
){}
void shade(const mgTLData& tld){}
void draw_STL(
	 const MGStl& stl, // 描画するMGStlオブジェクト
	 bool invokeNormal // glNormalを呼ぶかどうか
){}
void draw_points(
	const MGColor& boundary_color,
	const MGColor& inner_color,
	const std::vector<MGPosition>& ipos
){}
}//End of namespace mgGDL.

size_t MGObject::make_display_list(
	double span_length,//span length to approximate by polyline.
	int line_density//line density to draw surface in wire mode.
)const{	return dlist_name();}
void MGObject::delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be apppended.
)const{}

///Make a display list without color of this gel.
///Return is the display list name.
size_t MGObject::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{	return dlist_name()+1;}

size_t MGGroup::make_display_list(
	double span_length,//span length to approximate by polyline.
	int line_density//line density to draw surface in wire mode.
)const{	return dlist_name();}
///Make a display list without color of this gel.
///Return is the display list name.
size_t MGGroup::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{	return dlist_name()+1;}

void MGGroup::delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be apppended.
)const{}
void MGPlane::display_arrows()const{}
void MGStl::shade(
	double span_length	//Line segment span length.
)const{}
void MGStl::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGStl::display_arrows()const{}
void MGSurface::shade(
	double span_length	//Line segment span length.
)const{}
void MGFSurface::drawWireFS(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGShell::delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be apppended.
)const{}
void MGCurve::display_arrows()const{;}
void MGCurve::display_break_points()const{;};
void MGCurve::display_curvatures(
	double	scale,	///<scaling of the graph.
	int		density,///<densitiy of the graph.
	bool	use_radius///<true:radius display, false:curvature display.
)const{;};
void MGLBRep::display_control_polygon()const{}
void MGRLBRep::display_control_polygon()const{}
void MGSurface::display_arrows()const{}
void MGSBRep::display_control_polygon()const{}
void MGRSBRep::display_control_polygon()const{}
void MGFace::display_arrows()const{;}
void MGFace::display_control_polygon()const{}
void MGShell::display_arrows()const{}
void MGShell::display_control_polygon()const{}
MGGLAttrib* MGNullGLAttrib(long TID){return 0;}
void MGComplex::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGComplex::drawWire_in_star(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}

//Draw 3D point(vertex) in world coordinates.
//The object is converted to point(s) and is drawn.
void MGComplex::draw3DVertex(
)const{}
void MGComplex::draw3DVertex_in_star(
)const{}
void MGCellBase::drawWire_in_star(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGCellBase::draw3DVertex_in_star(
)const{}
void MGPVertex::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGPVertex::draw3DVertex(
)const{}

void MGBVertex::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGBVertex::draw3DVertex()const{}
void MGEdge::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{}
void MGEdge::draw3DVertex()const{}
void MGFace::draw3DVertex()const{}
void MGFace::shade(
	double span_length	//Line segment span length.
)const{}

//Shade the object in world coordinates.
void MGShell::shade(
	double span_length	//Line segment span length.
)const{}
void MGDraw_in_parameter_space(
	const MGObject& obj,	//The object to draw.
	double span_length		//Line segment span length.
){}

///Make a display list without color of this gel.
///Return is the display list name.
size_t MGFSurface::make_display_list_to_hilightFS(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{	return object_pointer()->dlist_name()+1;}

void MGShell::make_only_call_list(
	int mode	//=0 when wire, =1 when wire and no color, =2 when SHADING.
)const{}
size_t MGShell::make_display_list(
	double span_length,//span length to approximate by polyline.
	int line_density//line density to draw surface in wire mode.
)const{	return dlist_name();}
///Make a display list without color of this gel.
///Return is the display list name.
size_t MGShell::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{	return dlist_name()+1;}

void set_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit){mask|=bit;}
void reset_Amask(unsigned int& mask, MGGLAttrib::ATTRIB_MASK bit){mask&=~bit;}
MGColors::MGColors(){
	MGColors::m_colors=new MGColor*[MGColor::endID+1];
	MGColors::m_colors[0]=new MGColor(1.,1.,1.);
	for(size_t i=1; i<=MGColor::endID; i++)
		MGColors::m_colors[i]=new MGColor(1.,1.,1.);
}
MGColors::~MGColors(){
	for(int i=0; i<=MGColor::endID; i++)delete m_colors[i];
	delete[] m_colors;
}
namespace Gdiplus{
	int GdiplusStartup( unsigned long*, GdiplusStartupInput*, unsigned* ){return 0;}
	int GdiplusShutdown( unsigned long){return 0;}
};
void MGGroup::make_only_call_list_sub(
	size_t name,	//the name to glPushName for pick
	int mode,	//=0 when wire, =1 when wire and no color, =2 when SHADING.
	const std::vector<const MGGel*>* gels_to_delete
)const{}

#endif // defined MGCL_NO_MFC
