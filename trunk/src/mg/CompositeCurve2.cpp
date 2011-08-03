/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/CCisect_list.h"
#include "mg/CompositeCurve.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCompositeCurve Class.
//MGCompositeCurve is a composite of other leaf curves.
//Assumedly they are connected as C0 continuity. However, MGCompositeCurve
//does not check their continuity, but only put first or last as the user says
// (in connect_to_end or in connect_to_start), except two curves connecting
//are both MGLBRep or MGRLBRep. When the two connecting curves are both
//MGLBRep or MGRLBRep and they have more than C0 continuity,
//they are changed to one MGLBRep orMGRLBRep representation.
//Parameter ranges of the member curves are always continuous, from param_s() of the
//1st curve to param_e() of the last form MGCompositeCurve's paramter range.

//関数名：
//		common
//目的：
//		与えられた曲線と自身の交点もしくは共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			crv2,		(I/ )	与えられる曲線
//		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
//		MGCCisect_list&			isect		( /O)	交点
//				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
//							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
//戻り値：
//		3:交点も共通部分も求まった
//		2:交点のみが求まった
//		1:共通部分のみが求まった
//		0:交点も共通部分もなかった
//		-1:共通エッジの収束計算エラー
//		-2:共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()、をパラメータ範囲の収束計算の
//	誤差には、パラメータ範囲*rc_zero()を使用した
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan,
	MGCCisect_list& isect
)const{
	int obtainedIsect=0, obtainedCommon=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		MGCCisect_list isect2;
		int obtained2=(**i).common(crv2,vecComSpan2,isect2);
		if(obtained2<0)
			return obtained2;
		int ncom=vecComSpan2.size();
		if(ncom){
			obtainedCommon=1;
			for(int j=0; j<ncom; j++)
				vecComSpan.push_back(vecComSpan2[j]);
		}
		int nisect=isect2.size();
		if(nisect){
			obtainedIsect=2;
			isect.append(isect2);
		}
	}
	return obtainedCommon+obtainedIsect;
}

//関数名：
//		common
//目的：
//		与えられた曲線と自身の共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			crv2,		(I/ )	与えられる曲線
//		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
//				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
//							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
//戻り値：
//		共通部分の数:	共通部分が求まった
//		0:				共通部分がなかった
//		-1:				共通エッジの収束計算エラー
//		-2:				共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
//  パラメータ範囲*rc_zero()を使用した
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		int obtained2=(**i).common(crv2,vecComSpan2);
		if(obtained2<0)
			return obtained2;
		int ncom=vecComSpan2.size();
		for(int j=0; j<ncom; j++)
			vecComSpan.push_back(vecComSpan2[j]);

	}
	return vecComSpan.size()/4;
}

//一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
//は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線リストが返却される。
//costant offset curve. if the norm_vector is given, the positive offset direction decide
//to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGPvector<MGCurve> MGCompositeCurve::offset(
	double ofs_value,			//オフセット量
	const MGVector& norm_vector	//法線ベクトル
)const{	
return MGPvector<MGCurve>(static_cast<MGCurve*>(0));
	}

//可変オフセット関数
//オフセット量は空間次元1の線B表現で与えられる。
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
//は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線リストが返却される。
//valuable offset curve. if the norm_vector is given, the positive offset direction decide
//to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGPvector<MGCurve> MGCompositeCurve::offset(
	const MGLBRep& ofs_value_lb,					//空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector	//法線ベクトル
)const{
	return MGPvector<MGCurve>(static_cast<MGCurve*>(0));
}

//C2連続曲線の一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
//は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線が返却される。
//costant offset curve of C2 continuous curve. if the norm_vector is given, the positive offset direction
//decide to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGLBRep MGCompositeCurve::offset_c2(
	double ofs_value,			//オフセット量
	const MGVector& norm_vector	//法線ベクトル
)const{
	return MGLBRep();
}

//C2連続曲線の可変オフセット関数
//オフセット量は空間次元1の線B表現で与えられる。
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。ただし、曲率中心へ曲率半径以上のオフセット
//は行わない。トレランスはline_zero()を使用している。戻り値は、オフセット曲線が返却される。
//valuable offset curveof C2 continuous curve. if the norm_vector is given, the positive offset direction
//decide to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGLBRep MGCompositeCurve::offset_c2(
	const MGLBRep& ofs_value_lb,					//空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector	//法線ベクトル
)const{
	return MGLBRep();
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
MGSurface* MGCompositeCurve::sweep(
	const MGUnit_vector& uvec,			//Sweep Direction.
	double start_dist,					//distance to start edge.
	double end_dist			//distance to end edge.
)const{
	return new MGPlane();
}

//Obtain the projected curve of this onto the surface srf.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int MGCompositeCurve::project_onto_surface(
	const MGFSurface& srf,
	MGPvector<MGCurve>& vec_crv_uv,	
		//Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		srf.project(**i,vec_crv_uv,vec_crv,vec);
	return vec_crv.size();
}
