/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Default.h"
#include "mg/Ifstream.h"
#include "mg/Geometry.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/Point.h"
#include "mg/CompositeCurve.h"
#include "mg/MGStl.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//identify_type implements.
long MGPoint::identify_type()const{return MGPOINT_TID;}
long MGStraight::identify_type()const{return MGSTRAIGHT_TID;}
long MGEllipse::identify_type()const{return MGELLIPSE_TID;}
long MGLBRep::identify_type()const{return MGLBREP_TID;}
long MGRLBRep::identify_type()const{return MGRLBREP_TID;}
long MGSurfCurve::identify_type()const{return MGSRFCRV_TID;}
long MGTrimmedCurve::identify_type()const{return MGTRMCRV_TID;}
long MGCompositeCurve::identify_type()const{return MGCOMPCRV_TID;}
long MGPlane::identify_type()const{return MGPLANE_TID;}
long MGSphere::identify_type()const{return MGSPHERE_TID;}
long MGCylinder::identify_type()const{return MGCYLINDER_TID;}
long MGSBRep::identify_type()const{return MGSBREP_TID;}
long MGRSBRep::identify_type()const{return MGRSBREP_TID;}
long MGStl::identify_type()const{return MGSTL_TID;}

//------------------------------------------------
// MGGeometryのシリアライズ関数群

//メンバデータを書き込む関数
//MGGeometryのデータメンバが追加された場合にはここに処理を追加すること。
void MGGeometry::WriteMembers(MGOfstream& buf)const{
	MGObject::WriteMembers(buf);
	if(m_box) m_box->dump(buf);
	else MGDefault::null_box().dump(buf);
}

//メンバデータを読み出す関数
//MGGeometryのデータメンバが追加された場合にはここに処理を追加すること。
void MGGeometry::ReadMembers(MGIfstream& buf){
	MGObject::ReadMembers(buf);
	m_box=new MGBox();
	m_box->restore(buf);
	if(m_box->is_null()){
		delete m_box; m_box=0;
	}
}

//------------------------------------------------
// MGPointのシリアライズ関数群

//メンバデータを書き込む関数
void MGPoint::WriteMembers(MGOfstream& buf)const{
	MGGeometry::WriteMembers(buf);
	m_point.dump(buf);
}

//メンバデータを読み出す関数
void MGPoint::ReadMembers(MGIfstream& buf){
	MGGeometry::ReadMembers(buf);
	m_point.restore(buf);
}

//------------------------------------------------
// MGCurveのシリアライズ関数群

//メンバデータを書き込む関数
void MGCurve::WriteMembers(MGOfstream& buf)const{
	MGGeometry::WriteMembers(buf);
}

//メンバデータを読み出す関数
void MGCurve::ReadMembers(MGIfstream& buf){
	MGGeometry::ReadMembers(buf);
}

//------------------------------------------------
// MGStraightのシリアライズ関数群

//メンバデータを書き込む関数
void MGStraight::WriteMembers(MGOfstream& buf)const{
	MGCurve::WriteMembers(buf);
	m_root_point.dump(buf);
	m_direction.dump(buf);
	m_sparam.dump(buf);
	m_endparam.dump(buf);
}

//メンバデータを読み出す関数
void MGStraight::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_root_point.restore(buf);
	m_direction.restore(buf);
	m_sparam.restore(buf);
	m_endparam.restore(buf);
}

//------------------------------------------------
// MGEllipseのシリアライズ関数群

//メンバデータを書き込む関数
void MGEllipse::WriteMembers(MGOfstream& buf)const{
	MGCurve::WriteMembers(buf);
	m_center.dump(buf);
	m_normal.dump(buf);
	m_m.dump(buf);
	m_n.dump(buf);
	buf << m_r;
	buf<<m_prange[0];
	buf<<m_prange[1];
	buf<< m_circle;
	if(m_gprange){
		buf<<0x00000003L; //Flag of m_gprange data exists.
		buf<<m_gprange[0]<<m_gprange[1]<<m_gprange[2];
	}else{
		buf<<0x00000000L; //Flag of m_gprange does not data exists.
	}
}

//メンバデータを読み出す関数
void MGEllipse::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_center.restore(buf);
	m_normal.restore(buf);
	m_m.restore(buf);
	m_n.restore(buf);
	buf >> m_r;
	buf >>m_prange[0];
	buf >>m_prange[1];
	buf >> m_circle;
	long gprange;
	buf>>gprange;
	if(gprange){
		m_gprange=new double[3];
		buf>>m_gprange[0]>>m_gprange[1]>>m_gprange[2];
	}
}

//------------------------------------------------
// MGLBRepのシリアライズ関数群

//メンバデータを書き込む関数
void MGLBRep::WriteMembers(MGOfstream& buf)const{
	MGCurve::WriteMembers(buf);
	m_knot_vector.dump(buf);
	m_line_bcoef.dump(buf);
}

//メンバデータを読み出す関数
void MGLBRep::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_knot_vector.restore(buf);
	m_line_bcoef.restore(buf);
}

//------------------------------------------------
// MGRLBRepのシリアライズ関数群

//メンバデータを書き込む関数
void MGRLBRep::WriteMembers(MGOfstream& buf)const{
	MGCurve::WriteMembers(buf);
	m_line.WriteMembers(buf);
}

//メンバデータを読み出す関数
void MGRLBRep::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_line.ReadMembers(buf);
}

//------------------------------------------------
// MGSurfCurve のシリアライズ関数群

// メンバデータを書き込む関数
void MGSurfCurve::WriteMembers(MGOfstream& buf)const{
	MGCurve::WriteMembers(buf);
	buf.WritePointer(m_surface);
	m_curve.WriteMembers(buf);
}

void MGSurfCurve::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_surface=static_cast<MGSurface*>(buf.ReadPointer());
	m_curve.ReadMembers(buf);
}

//------------------------------------------------
// MGTrimmedCurve のシリアライズ関数群

// メンバデータを書き込む関数
void MGTrimmedCurve::WriteMembers(MGOfstream& buf)const{
	MGCurve::WriteMembers(buf);
	buf.WritePointer(m_curve);// pointer-basedなメンバデータの書き出し。
	m_range.dump(buf);
	if(m_sameRange) buf<<0x00000001L;
	else            buf<<0x00000000L;
}
void MGTrimmedCurve::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_curve=static_cast<MGCurve*>(buf.ReadPointer());
	m_range.restore(buf);
	long sameRange;
	buf>>sameRange;
	if(sameRange) m_sameRange=true;
	else          m_sameRange=false;
}

//------------------------------------------------
// MGCompositeCurve のシリアライズ関数群

//メンバデータを書き込む関数
void MGCompositeCurve::WriteMembers(MGOfstream& buf) const{
	size_t n=number_of_curves();
	buf<<n;
	for(size_t i=0; i<n; i++) buf.WritePointer(m_composite[i]);
}

//メンバデータを書き込む関数
void MGCompositeCurve::ReadMembers(MGIfstream& buf){
	size_t n;
	buf>>n;	//Number of element curves.
	MGCurve* crv;
	for(size_t i=0; i<n; i++){
		crv=static_cast<MGCurve*>(buf.ReadPointer());
		connect_to_end(crv);
	}
}

//------------------------------------------------
// MGSurface のシリアライズ関数群

//メンバデータを書き込む関数
//MGGeometryのデータメンバが追加された場合にはここに処理を追加すること。
void MGSurface::WriteMembers(MGOfstream& buf)const{
	MGGeometry::WriteMembers(buf);
}

//メンバデータを読み出す関数
//MGGeometryのデータメンバが追加された場合にはここに処理を追加すること。
void MGSurface::ReadMembers(MGIfstream& buf){
	MGGeometry::ReadMembers(buf);
}

//------------------------------------------------
// MGPlaneのシリアライズ関数群

//メンバデータを書き込む関数
void MGPlane::WriteMembers(MGOfstream& buf)const{
	MGSurface::WriteMembers(buf);
	m_normal.dump(buf);
	buf << m_d;
	m_root_point.dump(buf);
	m_uderiv.dump(buf);
	m_vderiv.dump(buf);
}

//メンバデータを読み出す関数
void MGPlane::ReadMembers(MGIfstream& buf){
	MGSurface::ReadMembers(buf);
	m_normal.restore(buf);
	buf >> m_d;
	m_root_point.restore(buf);
	m_uderiv.restore(buf);
	m_vderiv.restore(buf);
}
//------------------------------------------------
// MGSBRepのシリアライズ関数群

//メンバデータを書き込む関数
void MGSBRep::WriteMembers(MGOfstream& buf)const{
	MGSurface::WriteMembers(buf);
	m_surface_bcoef.dump(buf);
	m_uknot.dump(buf);
	m_vknot.dump(buf);
}

//メンバデータを読み出す関数
void MGSBRep::ReadMembers(MGIfstream& buf){
	MGSurface::ReadMembers(buf);
	m_surface_bcoef.restore(buf);
	m_uknot.restore(buf);
	m_vknot.restore(buf);
}

//------------------------------------------------
// MGRSBRepのシリアライズ関数群

//メンバデータを書き込む関数
void MGRSBRep::WriteMembers(MGOfstream& buf)const{
	MGSurface::WriteMembers(buf);
	m_surface.WriteMembers(buf);
}

//メンバデータを読み出す関数
void MGRSBRep::ReadMembers(MGIfstream& buf){
	MGSurface::ReadMembers(buf);
	m_surface.ReadMembers(buf);
}
