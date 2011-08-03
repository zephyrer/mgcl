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
 * �l�f�b�k�� B-Spline�Ȑ��AB-Spline�Ȗʂ���舵�����߂̂b�{�{��
 * �N���X���C�u�����ł��B�P�̂�NURBS�ȖʁA�g�������ꂽNURBS�ȖʁA����т�����
 * Winged Edge Data Structure�I�ȍ\���ɂ�茋�������V�F���\���܂ł̕\���\�͂�L���Ă��܂��B
 * MGCL�́A�O���t�B�b�N�ȕ\�������邽�߂ɂ� MFC(Microsoft Foundation Class) �� OpenGL ���A
 * �e�N�X�`���[�}�b�s���O�𗘗p���邽�߂ɂ�GDI-plus���A�����ă}�E�X�̃C�x���g�Ȃǂ�
 * �}���}�V���C���^�t�F�[�X�̂��߂�MFC�𗘗p���Ă��܂����A
 * �����A�����̋@�\�𗘗p���Ȃ��̂ł���΁A������r�����ĒP�Ȃ郉�C�u�����Ƃ��Ă����p�\�ŁA
 * �����ւ�Ɨ����̍������C�u�����ł�.
 * 
 * @section MGCL �l�f�b�k�̃N���X�A�֐��͉��L�̂悤�ȃ��W���[���ɕ��ނ���Ă��܂�
 * (1) Base Class
 * �@�@�l�f�b�k�̃x�[�X�ƂȂ�N���X�Q�ŁA�s�񉉎Z�A�x�N�^�A�e��R���e�i�A�{�b�N�X�g�Ȃǂ����S�ƂȂ�܂�.
 * 
 * (2) (Template) Functions or classe
 * �@�@�l�f�b�k�ŗ��p���Ă��鐔�l�v�Z�֐��ȂǗL�p�Ȋ֐����܂Ƃ߂Ă���܂�.
 * 
 * (3) Geometry(sub) classes
 * �@�@�_(MGPoint)�A���iMGCurve)�A��(MGSurface)�̊􉽕\���̂��߂̃N���X�� MGGeometry ��
 *     �T�u�N���X�ɂȂ�܂�.
 * 
 * (4) Topology(sub) classes
 * �@�@���_(MGPVertex, MGBVertex)�A�G�b�W(MGEdge)�A���[�v(MGLoop)�A�g�����Ȗ�(MGFace)�A
 *      �V�F��(MGShell)�̈ʑ��\���\���̂��߂̃N���X�� MGTopology �̃T�u�N���X�ɂȂ�܂�.
 * 
 * (5) GeoRelated classes
 * �@�@Geometery�N���X�̓��o�̓p�����[�^�ɗ��p�����e��N���X(MGCoons, MGPPRep, �[�������A
 *      �Ȃ�)�����ނ���܂�.
 * 
 * (6) TopoRelated classes
 * �@�@Topology�N���X�̓��o�̓p�����[�^�ɗ��p�����e��N���X(MGFOuterCurve, MGLPoint)��A
 *      �Ȗʂ̃g���������̃��[�e�B���e�B�N���X�� MGTrimLoop �Ȃǂ̃N���X�����ނ���܂�.
 * 
 * (7) Object Related classes
 * �@�@Geometry, Topology�ɋ���(MGFSurface)��A���̕��ނɓ��Ă͂܂�Ȃ��N���X(MGSTL)�A 
 *     MGObject �̃s�b�N�̂��߂̏��R���e�i�N���X�����ނ���܂�.
 * 
 * (8) Gel Related classes
 * �@�@ MGGel ��Group Element(Gel)�ƂȂ肤��N���X�̒��ۃN���X�ŁA���ׂĂ� MGObject �֘A�̍ō��ʂ�
 *     ���ۃN���X�ƂȂ�܂��B�O���[�v����Gel�̈ʒu��\������ MGGelPosition �Ȃǂ̃R���e�i�N���X�A
 *     MGGel �z���̂��ׂẴN���X���ނ̂��߂� MGAbstractGel, MGAbstractGels �Ȃǂ����ނ���܂�.
 * 
 * (9) File InputOutput classes
 * �@�@MGGel �z���̃N���X�͂܂��AMGCL�̗p�ӂ���t�@�C�����o�͂̑Ώۂł�����܂��B���̂��߂ɃN���X
 *     MGOfstream, MGIfstream, ������IGES�f�[�^�̓��o�͂̃N���X(MGIgesFstream, MGIgesIfstream,
 *     MGIgesOstream)�����ނ���܂�.
 * 
 * (10) Intersection Container classes
 * �@�@�_(MGPoint)�A��(MGCurve)�A��(MGSurface)�̊􉽕\�����m�A�܂��� MGEdge, MGLoop, MGFace, MGShell
 *     �Ȃǂ� Topology �N���X�Ƃ̌�_�Ƃ��̔z���\�����܂�.
 * 
 * (11) UseTessellation classes
 * �@�@��(MGSurface, MGFace)��\�����邽�߂ɂ͂������ׂ��ȎO�p�`�ŋߎ�����K�v������܂��B
 *     UseTessellation�͎O�p�`�ߎ��̂��߂̃c�[���N���X��񋟂��܂��B mgTLData, mgTLDataVector ��
 *     ����ŁA���̔z���ɁA���̋@�\�̎����̂��߂ɑ����̃N���X������܂����A��ʗ��p�҂͂��̂ӂ���
 *     �N���X�ŏ\���ł�.
 * 
 * (12) Display Handling classes
 * �@�@MGCL�ł͐}�`�̕\����OpenGL, �C���[�W�f�[�^�̎����̂��߂�GDI-plus�A�}���}�V���C���^�t�F�[�X��
 *     MFC(Micrsoft Foundation Class)�𗘗p���Ă��܂��BDisplay Handling Classes�ɂ͂l�f�b�k��
 *     �I�u�W�F�N�g�̕\���̂��߂̃N���X�ƕ\���̂��߂̃A�g���r���[�g�Ȃǂ����ނ���Ă��܂�.
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
// ���Ƃ���4.10�̏ꍇ�� MGCL0410�ɂȂ�B

///Get the MGCL_Version number.
MGEXTERN const char* MGCL_Version();

///Get the MGCL File validity.
MGEXTERN const char* MGCL_File_validity();

///  �� �l�̐ݒ�
const double mgPAI = 3.1415926535897932384626433833;

/// @var mgHALFPAI
/// @brief ��/2 �l�̐ݒ�
const double mgHALFPAI = mgPAI/2.;

/// 2.0 * �� �l�̐ݒ�
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
     MGINTERVAL_EMPTY,          ///<Empty interval. Interval �͋�W��
     MGINTERVAL_FINITE,         ///<Finite interval(above and below)
								///<����A�����Ƃ��L��
     MGINTERVAL_FINITE_ABOVE,   ///<Finite above and infinite below.
								///<����L���A��������
     MGINTERVAL_FINITE_BELOW,   ///<Finite below and infinite above.
								///<��������A�����L��
     MGINTERVAL_INFINITE        ///<Infinite below and above.
								///<����A�����Ƃ��ɖ��� 
};

///
///Relation type of plane and stright line.
///(Plane vs Plane, Plane vs S.Line, S.Line vs S.Line)
enum MGPSRELATION {
     MGPSREL_UNKNOWN,           ///<Unknown. �s��
     MGPSREL_TORSION,           ///<Torsion. �˂���
     MGPSREL_ISECT,				///<Intersection. ����
     MGPSREL_PARALLEL,          ///<Parallel. ���s
     MGPSREL_COIN,              ///<Coincidence. ���
     MGPSREL_VIRTUAL_ISECT      ///<Virtually intersection
								///<(Intersection at extended line).
								///<�w�蒼���̉�����̓_�ł̌���
};

///Curve type(�Ȑ��̎��)
enum MGCURVE_TYPE {
	MGCURVE_UNKNOWN,		///<Unkown. �s��
	MGCURVE_STRAIGHT,		///<MGStraight(Straight line). ����
	MGCURVE_ELLIPSE,		///<MGEllipse(Ellipse). �ȉ~
	MGCURVE_SPLINE,			///<B-Spline(MGLBRep). �X�v���C��
	MGCURVE_RSPLINE,		///<B-Spline(MGRLBRep)
	MGCURVE_SURFACE,		///<MGSurfCurve(Parameter line of a surface)
	MGCURVE_TRIMMED,		///<MGTrimmedCurve(Trimmed curve)
	MGCURVE_COMPOSITE,		///<MGCompositeCurve(Composite Curve)
	MGCURVE_USER1,			///<Auxiliary curve type 1
	MGCURVE_USER2,			///<Auxiliary curve type 2
	MGCURVE_BSUM			///<Boolean sum curve
};

///Ellipse type(�ȉ~�̎��)
enum MGELLIPSE_TYPE {
	 MGELLIPSE_EMPTY		///<Empty ellipse. ��
	,MGELLIPSE_SEGMENT		///<Segment ellipse. �Z�O�����g
	,MGELLIPSE_CLOSED		///<Closed(Whole) ellipse.�����ȉ~
};

///Straight line type(�����̎��)
enum MGSTRAIGHT_TYPE {
	 MGSTRAIGHT_EMPTY 		///<Empty. ��
	,MGSTRAIGHT_SEGMENT		///<Line segment. ����
	,MGSTRAIGHT_HALF_LIMIT	///<Half unlimit. ������
	,MGSTRAIGHT_UNLIMIT		///<Unlimit line for both direction. ��������
};

///Surface type(�Ȗʂ̎��)
enum MGSURFACE_TYPE {
	MGSURFACE_UNKNOWN,		///<Unknown. �s��
	MGSURFACE_PLANE,		///<Plane. ����
	MGSURFACE_CONE,			///<Cone. �~��
	MGSURFACE_SPHERE,		///<Sphere. ����
	MGSURFACE_TORUS,		///<Torus. �~��
	MGSURFACE_SPLINE,		///<Free form surface
							///<(Tensor product surface of LBRep) ���R�Ȗ�
	MGSURFACE_RSPLINE,		///<Free form surface
							///<(Rational B-Spline Surface)
	MGSURFACE_CYLINDER,		///<Cylinder surface(A special case of MGSURFACE_CONE)
	MGSURFACE_USER1,		///<Auxiliary surface type 1
	MGSURFACE_USER2,		///<Auxiliary surface type 2
	MGSURFACE_BSUM			///<Boolean sum surface
};

///Relation of curve and curve(�Ȑ��ƋȐ��̌�_�̊֌W)
enum MGCCRELATION {
	MGCCREL_UNKNOWN,	///<Unknown.
						///<���m�i���ׂĂ��Ȃ��āA�킩��Ȃ��A�ȉ��̂����ꂩ�j
	MGCCREL_ISECT,		///<Intersection. �����i�����A�ڂ���A��v�ȊO�j
	MGCCREL_NORMAL,		///<Intersection at right angle. ����
	MGCCREL_TANGENT,	///<Two curves are tangent. �ڂ���
	MGCCREL_COIN		///<Two curves are coincident. ��v����
};

///Relation of curve and surface(�Ȑ��ƋȖʂ̌�_�̊֌W)
enum MGCSRELATION {
	MGCSREL_UNKNOWN,///<Unknown. ���m
	MGCSREL_IN,		///<Intersection from inner of surface. �Ȗʂ̓�������̌�_
	MGCSREL_OUT,	///<ntersection from outer of surface. �Ȗʂ̊O������̌�_
	MGCSREL_IN_TAN,	///<Tangent from inner of surface.�Ȗʂ̓�������ڂ��Ă���
	MGCSREL_OUT_TAN,///<Tangent from outer of surface.�Ȗʂ̊O������ڂ��Ă���
	MGCSREL_COIN	///<Curve is included in surface. �Ȑ����ȖʂɊ܂܂�Ă���
};

///Relation of Surface and Surface(Surface��Surface�̌���̊֌W)
enum MGSSRELATION {
	MGSSREL_UNKNOWN,	///<Unknown.
						///<���m�i���ׂĂ��Ȃ��āA�킩��Ȃ��A�ȉ��̂����ꂩ�j
	MGSSREL_ISECT,		///<Intersection. �����i�ڂ���A��v�ȊO�j
	MGSSREL_TANGENT,	///<Tangent. �ڂ���
	MGSSREL_COIN		///<Coincident. ��v����
};

///End condition to get spline by interpolation.
enum MGENDCOND {
	MGENDC_UNKNOWN=0,	///< Unknown(usually not used).���m
	MGENDC_1D  =1,		///< 1st deravative provided.
	MGENDC_2D  =2,		///< 2nd deravative provided.
	MGENDC_NO  =3,		///< no end cond(only positional data)
	MGENDC_12D =4		///< both 1st and 2nd deravatives provided.
};

///a set of triangl type(�R�p�`���_���X�g�̃^�C�v)
enum mgTESTRIANG {
	mgTESTRIANG_UNKNOWN,	///<Unknown
	mgTESTRIANG_FAN,		///<like a csTriFanSet
	mgTESTRIANG_STRIP		///<like a csTriStripSet
};

///Tessellation ��subdivide���ꂽ�l�p�`��status.
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

///Dbug function. �f�o�b�O�֐�
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
