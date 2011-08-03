/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/*
 * Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno.
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
//
//  MGCL.h
//  MGCL.h provides enum variables.
/*! @mainpage Maestro's Geometry Class Library
 *
 *
 * ＭＧＣＬは B-Spline曲線、B-Spline曲面を取り扱うためのＣ＋＋の
 * クラスライブラリです。単体のNURBS曲面、トリムされたNURBS曲面、およびそれらを
 * Winged Edge Data Structure的な構造により結合したシェル構造までの表現能力を有しています。
 * MGCLは、グラフィックな表示をするためには MFC(Microsoft Foundation Class) と OpenGL を、
 * テクスチャーマッピングを利用するためにはGDI-plusを、そしてマウスのイベントなどの
 * マンマシンインタフェースのためにMFCを利用していますが、
 * もし、これらの機能を利用しないのであれば、これらを排除して単なるライブラリとしても利用可能で、
 * たいへん独立性の高いライブラリです.
 * 
 * @section MGCL ＭＧＣＬのクラス、関数は下記のようなモジュールに分類されています
 * (1) Base Class
 * 　　ＭＧＣＬのベースとなるクラス群で、行列演算、ベクタ、各種コンテナ、ボックス枠などが中心となります.
 * 
 * (2) (Template) Functions or classe
 * 　　ＭＧＣＬで利用している数値計算関数など有用な関数をまとめてあります.
 * 
 * (3) Geometry(sub) classes
 * 　　点(MGPoint)、線（MGCurve)、面(MGSurface)の幾何表現のためのクラスで MGGeometry の
 *     サブクラスになります.
 * 
 * (4) Topology(sub) classes
 * 　　頂点(MGPVertex, MGBVertex)、エッジ(MGEdge)、ループ(MGLoop)、トリム曲面(MGFace)、
 *      シェル(MGShell)の位相構造表現のためのクラスで MGTopology のサブクラスになります.
 * 
 * (5) GeoRelated classes
 * 　　Geometeryクラスの入出力パラメータに利用される各種クラス(MGCoons, MGPPRep, 端末条件、
 *      など)が分類されます.
 * 
 * (6) TopoRelated classes
 * 　　Topologyクラスの入出力パラメータに利用される各種クラス(MGFOuterCurve, MGLPoint)や、
 *      曲面のトリム処理のユーティリティクラスの MGTrimLoop などのクラスが分類されます.
 * 
 * (7) Object Related classes
 * 　　Geometry, Topologyに共通(MGFSurface)や、その分類に当てはまらないクラス(MGSTL)、 
 *     MGObject のピックのための情報コンテナクラスが分類されます.
 * 
 * (8) Gel Related classes
 * 　　 MGGel はGroup Element(Gel)となりうるクラスの抽象クラスで、すべての MGObject 関連の最高位の
 *     抽象クラスとなります。グループ内のGelの位置を表現する MGGelPosition などのコンテナクラス、
 *     MGGel 配下のすべてのクラス分類のための MGAbstractGel, MGAbstractGels などが分類されます.
 * 
 * (9) File InputOutput classes
 * 　　MGGel 配下のクラスはまた、MGCLの用意するファイル入出力の対象でもあります。このためにクラス
 *     MGOfstream, MGIfstream, そしてIGESデータの入出力のクラス(MGIgesFstream, MGIgesIfstream,
 *     MGIgesOstream)が分類されます.
 * 
 * (10) Intersection Container classes
 * 　　点(MGPoint)、線(MGCurve)、面(MGSurface)の幾何表現同士、または MGEdge, MGLoop, MGFace, MGShell
 *     などの Topology クラスとの交点とその配列を表現します.
 * 
 * (11) UseTessellation classes
 * 　　面(MGSurface, MGFace)を表示するためにはそれらを細かな三角形で近似する必要があります。
 *     UseTessellationは三角形近似のためのツールクラスを提供します。 mgTLData, mgTLDataVector が
 *     それで、この配下に、この機能の実現のために多くのクラスがありますが、一般利用者はこのふたつの
 *     クラスで十分です.
 * 
 * (12) Display Handling classes
 * 　　MGCLでは図形の表示にOpenGL, イメージデータの実現のためにGDI-plus、マンマシンインタフェースに
 *     MFC(Micrsoft Foundation Class)を利用しています。Display Handling ClassesにはＭＧＣＬの
 *     オブジェクトの表示のためのクラスと表示のためのアトリビュートなどが分類されています.
 *
 */
#ifndef _MGCL_HH_
#define _MGCL_HH_

/** @defgroup BASE Base Class
 *  @{
 */

class MGFPoint;

//
// The following macros are used to enable export/import.
// MGDECL for global functions and global valiable values in declaration.
// MGCLASS for MGCL classes.
// MGEXTERN for global functions and global valiable values.
// MGTEMPLATE for MGCL classes with template.
#if defined(MGCL_IMPORTS) || defined(MGCL_EXPORTS)
#	define MGCL_DLL
#	if defined(MGCL_EXPORTS)
		// Build DLL.
#		define MGDECL       __declspec(dllexport)
#		define MGCLASS      __declspec(dllexport)
#		define MGEXTERN     extern __declspec(dllexport)
#		define MGTEMPLATE
#	else
		// Import DLL.
#		define MGDECL       __declspec(dllimport)
#		define MGCLASS      __declspec(dllimport)
#		define MGEXTERN     extern __declspec(dllimport)
#		define MGTEMPLATE   extern
#	endif
#else
	// Not DLL.
#	define MGDECL
#	define MGCLASS
#	define MGEXTERN         extern
#	define MGTEMPLATE

#endif	// MGCL_IMPORTS || MGCL_EXPORTS

#include <iosfwd>
#include <assert.h>
#include <vector>

MGEXTERN const char* _MGCL_VER;
MGEXTERN const char* _MGCL_FILE;
// Maestro's Geometry Classes Library version 5.30
// <<<<<<<<<<4/20/2001>>>>>>>>>>>> //
// たとえば4.10の場合は MGCL0410になる。

///Get the MGCL_Version number.
MGEXTERN const char* MGCL_Version();

///Get the MGCL File validity.
MGEXTERN const char* MGCL_File_validity();

///  π 値の設定
const double mgPAI = 3.1415926535897932384626433833;

/// @var mgHALFPAI
/// @brief π/2 値の設定
const double mgHALFPAI = mgPAI/2.;

/// 2.0 * π 値の設定
const double mgDBLPAI = mgPAI*2.;

///
///Infinite value Definition. This value is used in MGEReal class
///to identify infinite value of double.
///
const double mgInfiniteVal = 1.e+20;

///Infinite type
enum MGINFINITE_TYPE {
     MGINFINITE_MINUS=-1,      ///< Minus infinite
	 MGINFINITE_FINITE=0,      ///< Finite
     MGINFINITE_PLUS=1,		   ///< Plus infinite
};

///MGInterval type
enum MGINTERVAL_TYPE {
     MGINTERVAL_EMPTY,          ///<Empty interval. Interval は空集合
     MGINTERVAL_FINITE,         ///<Finite interval(above and below)
								///<上方、下方とも有限
     MGINTERVAL_FINITE_ABOVE,   ///<Finite above and infinite below.
								///<上方有限、下方無限
     MGINTERVAL_FINITE_BELOW,   ///<Finite below and infinite above.
								///<上方無限、下方有限
     MGINTERVAL_INFINITE        ///<Infinite below and above.
								///<上方、下方ともに無限 
};

///
///Relation type of plane and stright line.
///(Plane vs Plane, Plane vs S.Line, S.Line vs S.Line)
enum MGPSRELATION {
     MGPSREL_UNKNOWN,           ///<Unknown. 不明
     MGPSREL_TORSION,           ///<Torsion. ねじれ
     MGPSREL_ISECT,				///<Intersection. 交差
     MGPSREL_PARALLEL,          ///<Parallel. 平行
     MGPSREL_COIN,              ///<Coincidence. 包含
     MGPSREL_VIRTUAL_ISECT      ///<Virtually intersection
								///<(Intersection at extended line).
								///<指定直線の延長上の点での交差
};

///Curve type(曲線の種類)
enum MGCURVE_TYPE {
	MGCURVE_UNKNOWN,		///<Unkown. 不明
	MGCURVE_STRAIGHT,		///<MGStraight(Straight line). 直線
	MGCURVE_ELLIPSE,		///<MGEllipse(Ellipse). 楕円
	MGCURVE_SPLINE,			///<B-Spline(MGLBRep). スプライン
	MGCURVE_RSPLINE,		///<B-Spline(MGRLBRep)
	MGCURVE_SURFACE,		///<MGSurfCurve(Parameter line of a surface)
	MGCURVE_TRIMMED,		///<MGTrimmedCurve(Trimmed curve)
	MGCURVE_COMPOSITE,		///<MGCompositeCurve(Composite Curve)
	MGCURVE_USER1,			///<Auxiliary curve type 1
	MGCURVE_USER2,			///<Auxiliary curve type 2
	MGCURVE_BSUM			///<Boolean sum curve
};

///Ellipse type(楕円の種類)
enum MGELLIPSE_TYPE {
	 MGELLIPSE_EMPTY		///<Empty ellipse. 空
	,MGELLIPSE_SEGMENT		///<Segment ellipse. セグメント
	,MGELLIPSE_CLOSED		///<Closed(Whole) ellipse.閉じた楕円
};

///Straight line type(直線の種類)
enum MGSTRAIGHT_TYPE {
	 MGSTRAIGHT_EMPTY 		///<Empty. 空
	,MGSTRAIGHT_SEGMENT		///<Line segment. 線分
	,MGSTRAIGHT_HALF_LIMIT	///<Half unlimit. 半直線
	,MGSTRAIGHT_UNLIMIT		///<Unlimit line for both direction. 無限直線
};

///Surface type(曲面の種類)
enum MGSURFACE_TYPE {
	MGSURFACE_UNKNOWN,		///<Unknown. 不明
	MGSURFACE_PLANE,		///<Plane. 平面
	MGSURFACE_CONE,			///<Cone. 円錐
	MGSURFACE_SPHERE,		///<Sphere. 球面
	MGSURFACE_TORUS,		///<Torus. 円環面
	MGSURFACE_SPLINE,		///<Free form surface
							///<(Tensor product surface of LBRep) 自由曲面
	MGSURFACE_RSPLINE,		///<Free form surface
							///<(Rational B-Spline Surface)
	MGSURFACE_CYLINDER,		///<Cylinder surface(A special case of MGSURFACE_CONE)
	MGSURFACE_USER1,		///<Auxiliary surface type 1
	MGSURFACE_USER2,		///<Auxiliary surface type 2
	MGSURFACE_BSUM			///<Boolean sum surface
};

///Relation of curve and curve(曲線と曲線の交点の関係)
enum MGCCRELATION {
	MGCCREL_UNKNOWN,	///<Unknown.
						///<未知（調べていなくて、わからない、以下のいずれか）
	MGCCREL_ISECT,		///<Intersection. 交差（直交、接する、一致以外）
	MGCCREL_NORMAL,		///<Intersection at right angle. 直交
	MGCCREL_TANGENT,	///<Two curves are tangent. 接する
	MGCCREL_COIN		///<Two curves are coincident. 一致する
};

///Relation of curve and surface(曲線と曲面の交点の関係)
enum MGCSRELATION {
	MGCSREL_UNKNOWN,///<Unknown. 未知
	MGCSREL_IN,		///<Intersection from inner of surface. 曲面の内側からの交点
	MGCSREL_OUT,	///<ntersection from outer of surface. 曲面の外側からの交点
	MGCSREL_IN_TAN,	///<Tangent from inner of surface.曲面の内側から接している
	MGCSREL_OUT_TAN,///<Tangent from outer of surface.曲面の外側から接している
	MGCSREL_COIN	///<Curve is included in surface. 曲線が曲面に含まれている
};

///Relation of Surface and Surface(SurfaceとSurfaceの交線の関係)
enum MGSSRELATION {
	MGSSREL_UNKNOWN,	///<Unknown.
						///<未知（調べていなくて、わからない、以下のいずれか）
	MGSSREL_ISECT,		///<Intersection. 交差（接する、一致以外）
	MGSSREL_TANGENT,	///<Tangent. 接する
	MGSSREL_COIN		///<Coincident. 一致する
};

///End condition to get spline by interpolation.
enum MGENDCOND {
	MGENDC_UNKNOWN=0,	///< Unknown(usually not used).未知
	MGENDC_1D  =1,		///< 1st deravative provided.
	MGENDC_2D  =2,		///< 2nd deravative provided.
	MGENDC_NO  =3,		///< no end cond(only positional data)
	MGENDC_12D =4		///< both 1st and 2nd deravatives provided.
};

///a set of triangl type(３角形頂点リストのタイプ)
enum mgTESTRIANG {
	mgTESTRIANG_UNKNOWN,	///<Unknown
	mgTESTRIANG_FAN,		///<like a csTriFanSet
	mgTESTRIANG_STRIP		///<like a csTriStripSet
};

///Tessellation のsubdivideされた四角形のstatus.
enum MGRECT_STATUS{
	MGRECT_UNKNOWN,
	MGRECT_IN,		///<Whole rectangle is inside the Face.
	MGRECT_OUT,		///<Whole rectangle is outside the Face.
	MGRECT_ON,		///<Some trimming curve is crossing the rectangle.
					///<So part is inside, and part is outside the face.
	MGRECT_OVER,	///<In the rectangle whole inner boundary loop is included.
	MGRECT_ONANDOVER///<In the rectangle whole inner boundary loop is included, and
					///<Some trimming curve is crossing the rectangle.
};

///Dbug function. デバッグ関数
MGDECL std::ostream& operator<< (std::ostream& out, MGINTERVAL_TYPE rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGPSRELATION rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGCURVE_TYPE rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGELLIPSE_TYPE rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGSTRAIGHT_TYPE rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGSURFACE_TYPE rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGCCRELATION rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGCSRELATION rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGSSRELATION rel);
MGDECL std::ostream& operator<< (std::ostream& out, MGENDCOND rel);

///MGCL namespace defines varialbes without prefix mg or MG.
///From historical reasons, varialbes with prefix mg or MG are not included in MGCL namespace.
namespace MGCL{

///Tessellation parameter to select fan kind for the tessellation.
enum fan_kind{
	SINGLE_TRIANGLE, ///< 1 triangle/FAN(default) and STRIP for as many as posible triangles,
		///<STRIP triangles may cover multiple rectangles.

	MULTIPLE_TRIANGLES, ///< as many triangles as possible/FAN and STRIP for as many as posible triangles,
		///<STRIP triangles may cover multiple rectangles.

	SINGLE_TRIANGLE_NO_STRIP,///<SINGLE_TRIANGLE, but STRIP triangles cover only one tessellated rectagle.
	MULTIPLE_TRIANGLES_NO_STRIP,///<MULTIPLE_TRIANGLES, but STRIP triangles cover only one tessellated rectagle.
};

enum SURFACE_CURVATURE_KIND{
	GAUSSIAN_CURVATURE = 0, ///< Gaussian curvature
	MEAN_CURVATURE = 1,  ///< mean curvature
	MINIMUM_CURVATURE = 2,   ///< the minimum curvature
	MAXIMUM_CURVATURE = 3    ///< the maximum curvature
};

///
///Display mode of MGObject.
///SHADING and WIRE_AND_SHADING are valid only for MGObject of manifold dimension 2.
///
enum VIEWMODE{
	WIRE=0,		///< wire frame mode
	SHADING=1,	///< surface mode
	WIRE_AND_SHADING=2 ///< wire and surface mode
};

///convert the angel unit from degree to radian.
inline double degree_to_radian(double degree){
	const double coef = mgPAI / 180.;
	return degree * coef;
}
///convert the angel unit from radian to degree.
inline double radian_to_degree(double radian){
	const double coef = 180. / mgPAI;
	return radian * coef;
}

///Start up the MGCL.
///This is necessary only when MGIgesxxxx class or MGImage or MGTexturexxxxx class is
///to use. Before use of GDIplus, GdiStartUp is necessary, this start_up will do it.
///
void start_up(
	bool need_to_GdiStartUp=false	///<True if GdiplusStartUp is necessary.
);

///Shut down the MGCL.
void shut_down();

static unsigned long m_gdiplusToken;///<a token to pass GdiplusShutdown.
							///<Initialized at GdiplusStartup.
static bool m_gdiplus_initialized;///<Indicates if MGCL startup Gdiplus.

///Compute the difference of min and max of the three doubles a1, a2, and a3.
double Max3(double a1, double a2, double a3);

}

/** @} */ // end of BASE group
#endif
