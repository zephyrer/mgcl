/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGGLDrawFunc_HH_
#define _MGGLDrawFunc_HH_

#include <vector>
#include "MGCLStdAfx.h"

class MGBPointSeq;
class MGBox;
class MGVector;
class MGCurve;
class MGFSurface;
class MGPosition;
class MGSPointSeq;
class mgTLData;
class mgTLDataVector;
class MGTextureImage;
class MGTextureImage2;
class MGTextureImages;
class mgTLTriangle;
class MGStl;
class MGColor;

/** @addtogroup DisplayHandling
 *  @{
 */

///In the namespace mgGDL, usefule display functions are defined.
///Originally mgGDL is inteded not to depend on OpenGL.
///That is, if function body is exchanged to other implemetation than OpenGL,
///users of mgGDL are not affected.
namespace mgGDL{

MGDECL void MGGLBeginLINE_STRIP();

MGDECL void MGGLVertex(double x, double y, double z);

MGDECL void MGGLEnd();

///Draw an arrow symbol with implementation of OpenGL.
///pos[0] is the origin of the arrow, pos[1] is the top of the arrow,
///pos[2], [3] are two bottoms of arrowhead.
MGDECL void MGDrawArrow(MGPosition pos[4]);

///Draw arrow symbols of curve with implementation of OpenGL.
///ndiv is the number of arrows for the curve.
MGDECL void MGDrawArrow(const MGCurve& curve, int ndiv);

///Draw arrow symbols of curve with implementation of OpenGL.
///arrows of udiv*vdiv are displyed on the surface.
MGDECL void MGDrawArrow(const MGFSurface& surf, int udiv, int vdiv);

///Draw an object of class MGBox, by wireframe.
MGDECL void MGDrawBox(const MGBox&);

///Draw a control polygon.
///Lines can be drawn between point[i-1] and point[i], for i=1, ..., length()-1.
MGDECL void MGDrawPointSeq(
	const MGBPointSeq& bp,	///<input point sequence.
	bool draw_points=true,	///<True if points be drawn.
	bool dotted=true		///<true if dotted lines be drawn.
);

///Draw control polygon network.
///Lines can be drawn between sp(i,j,.).
MGDECL void MGDrawPointSeq(
	const MGSPointSeq&	sp,
	bool draw_points=true,	///<True if points be drawn.
	bool dotted=true		///<true if dotted lines be drawn.
);

///Draw curvature variation graph so-called Hige.
MGDECL void MGDrawCurvaGraph(const MGCurve& curve, double scale, int density, bool radius);

///Draw a point using openGL functions.
MGDECL void MGDrawPoint(double x,double y,double z);
MGDECL void MGDrawPoint(const MGPosition& pos);

///Draw a polyline using openGL functions.
///These three functions do not set the line attributes. So, if necessary,
///line attrubutes must be set before invoking the functions.
///When cloded=true, 1st and last points will be connected.
MGDECL void MGDrawPolyline(const MGBPointSeq& line, bool closed=false);
MGDECL void MGDrawPolyline(const std::vector<MGPosition>& line, bool closed=false);
MGDECL void MGDrawStraight(const MGPosition& end, const MGPosition& start);

///Generate the OpenGL display list for GL3ParamView's parameter rectangle.
MGDECL void param_rectangle(
	const MGBox& param_box	///<Parameter range.
);

///OpenGL shading display of a tesselated data tld.
MGDECL void shade(const mgTLData& tld);

///OpenGL display for the tessellation lines drawn in world view.
MGDECL void mgGDL_WTess(
	const mgTLData& tld	///<tessellation data.
);

///OpenGL declaration for shading.
///setup_easy_shade calls the follows:
///1. glEnable for GL_LIGHT0, GL_LIGHTING, GL_COLOR_MATERIAL.
///2. glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
///3. glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
///4. glDisable(GL_CULL_FACE);/// 両面を描画
///5. glColorMaterial(GL_FRONT_AND_BACK) for GL_AMBIENT & GL_DIFFUSE.
void setup_easy_shade(
);

///OpenGL declaration for shading. setup_wire_drawing calls the follows:
///	glDisable(GL_COLOR_MATERIAL);/// ライトをオフ
///	glDisable(GL_LIGHTING);/// ライティングオフ	
///	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);/// ワイヤーフレームを描画		
void setup_wire_drawing(
);

///OpenGL display for the tessellation lines drawn in parameter view.
MGDECL void mgGDL_PTess(
	const mgTLData& tld	///<tessellation data.
);

///Renders curvatures mapping that comes into colorful image.
///A point whose curvature is within [lower, upper], the color varies.
void draw_surface_curvature_data(
	const mgTLData& tld,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper, ///<minimum and maximum value of the curvatures of the kind.
		///<Whne lower>=upper, lower is set as the minimum value and upper is set
		///<as the maximum value out of all the curvatures.
	double* real_lower=0,	///<When not null, minimum curvature of the kind will be output.
	double* real_upper=0	///<When not null, maximum curvature of the kind will be output.
);

/// Renders curvatures mapping that comes into colorful image.
/// A point whose curvature is within [lower, upper], the color varies.
void draw_surface_curvature(
	const mgTLDataVector& tldvec,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper, ///<minimum and maximum value of the curvatures of the kind.
		///<Whne lower>=upper, lower is set as the minimum value and upper is set
		///<as the maximum value out of all the curvatures.
	double* real_lower=0,	///<When not null, minimum curvature of the kind will be output.
	double* real_upper=0	///<When not null, maximum curvature of the kind will be output.
);

///Set up the environment for texture draw.
void draw_texture_image_set_up(
	const double back_color[4]	///<back ground color data.
);

///End up process for the environment for texture draw.
void draw_texture_image_end_up();

/// MGStlオブジェクトを描画する関数
/// glBiginからglEndまでを実行する
void draw_STL(
	 const MGStl& stl,///< 描画するMGStlオブジェクト
	 bool invokeNormal=false///< glNormalを呼ぶかどうか
);

///draw poits sequence ipos with 2 colors, inner and outer.
void draw_points(
	const MGColor& boundary_color,
	const MGColor& inner_color,
	const std::vector<MGPosition>& ipos
);

}

/** @} */ // end of DisplayHandling group
#endif