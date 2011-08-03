/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// tangentCurve.cpp : MGCurveとMGSurfaceからTangent Planeを作る関数
//
#include "MGCLStdAfx.h"
#include "mg/LBRep.h"
#include "mg/Surface.h"
#include "mg/Position.h"
#include "mg/SurfCurve.h"
#include "mg/Interval.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGLBRep TP_from_world_curve(
	const MGSurface& srf,	//接続面
	const MGCurve& crv,		//接続曲線
	size_t order,			//作成するカーブのオーダ
	int& error)	//エラーコード  0:OK, !=0 error code when constructing LBRep.
{
	if(order < 2)return MGLBRep();
	const MGKnotVector& tempKnotVector = crv.knot_vector();
	MGKnotVector knotVector(tempKnotVector, order);	//指定オーダーのノットベクトルに作り替える
	size_t i;
	for(i = tempKnotVector.order() - 1; i < tempKnotVector.bdim(); i++){
		double	tpara = 0.0,					//テンポラリ
				spara = tempKnotVector(i),		//スパンの始点
				epara = tempKnotVector(i + 1);	//スパンの終点
		if(epara - spara < crv.param_error())continue;	//マルチノットのときの処理(RLBRepのみ)
		//1スパンの分割数を決定する
		MGInterval interval(spara, epara);
		int ndiv = crv.calc_div_num(interval);			//オフセット曲線の分割数
		double shortspan = (epara - spara) / ndiv;
		tpara = spara;
		for (int j = 0; j < ndiv; j++){knotVector.add_data(tpara); tpara += shortspan;}
	}

	// onを使うより perps_guessを使ったほうが早いので
	// 先に第一点を求めておく
	MGPosition guess_prm(2), on_prm(2);
	MGPosition crv_pt(crv.start_point()), srf_tan(srf.sdim());
	MGPosition uv_min(srf.param_s_u(), srf.param_s_v()), 
			uv_max(srf.param_e_u(), srf.param_e_v());
	srf.on(crv_pt, guess_prm);
	//制御点を生成する
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);
	int len = dataPoint.length();
	MGBPointSeq bp1(len, srf.sdim());
	for(int j = 0; j < len; j++){
		srf.perp_guess(uv_min, uv_max, 
			crv.eval(dataPoint(j)),guess_prm, on_prm);
		MGPosition pos = srf.unit_normal(on_prm);
		guess_prm = on_prm;
		bp1.store_at(j, pos);
	}

	//精度十分の曲線を生成する
	MGLBRep brep;
	brep = MGLBRep(dataPoint, bp1, knotVector, error);
	//余分なKnotを削除
	MGTolerance::push();
	// 2.0* sinθ/2 ≒ θ とする。
	MGTolerance::set_line_zero(MGTolerance::angle_zero());
	brep.remove_knot();
	MGTolerance::pop();
	return brep;
}

MGLBRep TP_from_parameter_curve(
	const MGSurface& srf,	//接続面
	const MGCurve& pcrv,	//接続曲線(接続面上のパラメータカーブ)
	const MGCurve& wcrv,	//接続曲線(世界座標上のカーブ)
	size_t order,			//作成するカーブのオーダ
	int& error)	//エラーコード  0:OK, !=0 error code when constructing LBRep.
{
	// WorldCurveから十分細かいKnotVectorを作成する。
	if(order < 2)return MGLBRep();
	const MGKnotVector& tempKnotVector = wcrv.knot_vector();
	MGKnotVector knotVector(tempKnotVector, order);	//指定オーダーのノットベクトルに作り替える
	size_t i;
	for(i = tempKnotVector.order() - 1; i < tempKnotVector.bdim(); i++){
		double	tpara = 0.0,					//テンポラリ
				spara = tempKnotVector(i),		//スパンの始点
				epara = tempKnotVector(i + 1);	//スパンの終点
		if(epara - spara < wcrv.param_error())continue;	//マルチノットのときの処理(RLBRepのみ)
		//1スパンの分割数を決定する
		MGInterval interval(spara, epara);
		int ndiv = wcrv.calc_div_num(interval);			//オフセット曲線の分割数を代用
		double shortspan = (epara - spara) / ndiv;
		tpara = spara;
		for (int j = 0; j < ndiv; j++){knotVector.add_data(tpara); tpara += shortspan;}
	}

	//制御点を生成する
	// SurfCurve(src,pcrv)のon()で
	// wcrvの或るパラメータtに対応するpcrvのパラメータ値ptを求め
	// ptからSurfaceのuv値およびNormalを求める
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);
	int len = dataPoint.length();
	MGBPointSeq bp1(len, srf.sdim());
	MGSurfCurve scrv(srf, pcrv);
	for(int j = 0; j < len; j++){
		double pt;
		scrv.on(wcrv.eval(dataPoint(j)),pt);
		MGPosition pos = srf.unit_normal(pcrv.eval(pt));
		bp1.store_at(j, pos);
	}

	// 冗長なLBRepを作る
	MGLBRep brep;
	brep = MGLBRep(dataPoint, bp1, knotVector, error);

	//余分なノットを落とす
	MGTolerance::push();
	// 2.0* sinθ/2 ≒ θ とする。
	MGTolerance::set_line_zero(MGTolerance::angle_zero());
	brep.remove_knot();
	MGTolerance::pop();
	return brep;
}
