/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Knot.h"
#include "mg/Curve.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Straight.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of curve offset.

#define POW_LOW 0.34	//最初にオフセットさせるポイント数を決める係数(速度重視)
#define POW_HIGH 0.50   //最初にオフセットさせるポイント数を決める係数(精度重視)
#define NUM_DIV 50      //1スパンの分割数がこれを越えたときPOW_LOWを使用するようにする

//向きが同じ2本のB表現曲線を接続する(同じ種類のとき)
MGLBRep* join2LBRep(const MGLBRep& crv1, const MGLBRep& crv2){
	int which, cont;
	double ratio;
	cont = crv1.continuity(crv2, which, ratio);
	if(cont < 0 || which == 0 || which == 3)
		return NULL;

	return new MGLBRep(crv1, cont, which, crv2);
}

//一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線リストが返却される。
MGPvector<MGCurve> MGCurve::offset(
	double ofs_value,			//オフセット量
	const MGVector& norm_vector	//法線ベクトル
)const{
	int error;
	MGBPointSeq bp1(2, 1);					//ofs_value一定の直線を生成する
	bp1.store_at(0, &ofs_value);
	bp1.store_at(1, &ofs_value);
	MGLBRep ofs_value_lb(bp1, error, 2);	//オーダー２の直線を作る
	if(error){
		return MGPvector<MGCurve>();
	}
	ofs_value_lb.change_range(param_s(), param_e());
	return offset(ofs_value_lb, norm_vector);	//可変オフセットを使用する
}

//可変オフセット関数
//オフセット量は空間次元1の線B表現で与えられる。
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線リストが返却される。
MGPvector<MGCurve> MGCurve::offset(
	const MGLBRep& ofs_value_lb,	//空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector		//法線ベクトル
)const{
	MGPvector<MGCurve> ofs_crvl;
	if(norm_vector == mgNULL_VEC){
		offset_proc(ofs_value_lb, ofs_crvl);
	}else{
		offset_norm_proc(ofs_value_lb, norm_vector, ofs_crvl);
	}
	return ofs_crvl;
}

//一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線が返却される。
MGPvector<MGCurve> MGStraight::offset(
	double ofs_value,			//オフセット量
	const MGVector& norm_vector	//法線ベクトル
)const{
	MGPvector<MGCurve> ofs_crvl;
	MGStraight *st1 = clone();
	if(norm_vector == mgNULL_VEC){
		*st1 += MGUnit_vector() * ofs_value;
	}else{
		MGUnit_vector tangent = st1->direction(param_s()), ofs_dir;
		ofs_dir = norm_vector * tangent;
		*st1 += ofs_dir * ofs_value;
	}
	ofs_crvl.push_back(st1);
	return ofs_crvl;
}

//法線ベクトルが指定されているときの処理
int MGCurve::offset_norm_proc(
	const MGLBRep& ofs_value_lb,	//オフセット量
	const MGVector& norm_vector,	//法線ベクトル
	MGPvector<MGCurve>& ofs_crvl	//オフセットカーブリスト
)const{
	//曲線を折れがあったら（multiple knot）で分割する
	MGPvector<MGCurve> crv_list, tmp_list;
	int nDiv = divide_multi_ofs(ofs_value_lb, crv_list);

	ofs_crvl.clear();
	for(int i = 0;i < nDiv; i++){
		MGLBRep tmp_brep;
		if(!crv_list[i]->offset_norm_c2_proc(ofs_value_lb, norm_vector, tmp_brep))return false;
		tmp_list.push_back(tmp_brep.clone());
	}
	int num = join(tmp_list, ofs_crvl);
	return num;
}

//法線ベクトルが指定されていないときの処理
int MGCurve::offset_proc(
	const MGLBRep& ofs_value_lb,		//オフセット量
	MGPvector<MGCurve>& ofs_crvl	//オフセットカーブリスト
)const{
	//曲線を折れ（multiple knot）で分割する
	MGPvector<MGCurve> crv_list, tmp_list;
	int nDiv = divide_multi_ofs(ofs_value_lb, crv_list);

	ofs_crvl.clear();
	MGUnit_vector T, B, preN;	//preN : 前回のノーマルベクトル
	int freverse = 0;			//向きを逆にしているというフラグ
	double curvature, torsion;
	Frenet_frame(param_s(), T, preN, B, curvature, torsion);	//カーブ始点のノーマルを最初に与える

	for(int i = 0;i < nDiv; i++){
		MGLBRep tmp_brep;
		if(!crv_list[i]->offset_c2_proc(ofs_value_lb, tmp_brep, preN, freverse))return false;
		tmp_list.push_back(tmp_brep.clone());
	}
	int num = join(tmp_list, ofs_crvl);
	return num;
}

//C2連続曲線の一定オフセット関数
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線が返却される。
MGLBRep MGCurve::offset_c2(
	double ofs_value,				//オフセット量
	const MGVector& norm_vector		//法線ベクトル
)const{
	int error;
	MGBPointSeq bp1(2, 1);					//ofs_value一定の直線を生成する
	bp1.store_at(0, &ofs_value);
	bp1.store_at(1, &ofs_value);
	MGLBRep ofs_value_lb(bp1, error, 2);	//オーダー２の直線を作る
	if(error)return MGLBRep();
	ofs_value_lb.change_range(param_s(), param_e());
	return offset_c2(ofs_value_lb, norm_vector);	//可変オフセットを使用する
}

//C2連続曲線の可変オフセット関数
//オフセット量は空間次元1の線B表現で与えられる。
//オフセット方向は、法線方向から見て入力曲線の進行方向左側を正とする。
//法線ベクトルがヌルの場合、始点において曲率中心方向を正とする。
//ただし、曲率中心へ曲率半径以上のオフセットは行わない。トレランスはline_zero()を使用している。
//戻り値は、オフセット曲線が返却される。
MGLBRep MGCurve::offset_c2(
	const MGLBRep& ofs_value_lb,	//空間次元１の線B表現で示したオフセット量
	const MGVector& norm_vector		//法線ベクトル
)const{
	int rc = 0;
	MGLBRep ofs_brep;
	if(norm_vector == mgNULL_VEC){
		int freverse = 0;			//向きを逆にしているというフラグ
		MGUnit_vector T, B, preN;	//preN : 前回のノーマルベクトル
		double curvature, torsion;
		Frenet_frame(param_s(), T, preN, B, curvature, torsion);	//カーブ始点のノーマルを最初に与える
		rc = offset_c2_proc(ofs_value_lb, ofs_brep, preN, freverse);
	}else{
		rc = offset_norm_c2_proc(ofs_value_lb, norm_vector, ofs_brep);
	}
	if(!rc)return MGLBRep();
	ofs_brep.remove_knot();
	return ofs_brep;
}

//法線ベクトルが指定されているC2連続曲線のオフセット
int MGCurve::offset_norm_c2_proc(
	const MGLBRep& ofs_value_lb,//オフセット量
	const MGVector& norm_vector,//法線ベクトル
	MGLBRep& ofs_brep			//オフセットカーブ
)const{
	MGKnotVector knotVector = offset_make_knotvector(ofs_value_lb);	//十分分割したノットベクトルを求める
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);

	//制御点を生成する
	int len = dataPoint.length();
	MGBPointSeq bp1(len, sdim());
	for(int i = 0; i < len; i++){
		MGUnit_vector T, N, B, ofs_dir;
		double curvature, torsion, ofs_value = ofs_value_lb.eval_position(dataPoint(i)).ref(0);
		Frenet_frame(dataPoint(i), T, N, B, curvature, torsion);
		ofs_dir = norm_vector * T;
		double cosine = N % ofs_dir;
		int fneg = (cosine < 0);	//コサインの正負フラグ
		if(!MGMZero(curvature))		//ゼロ割チェック
			if((ofs_value * (-2. * fneg + 1)) > ((1. / curvature) / fabs(cosine)))return false;
		MGPosition pos;
		pos = eval(dataPoint(i));
		pos += ofs_dir * ofs_value;
		bp1.store_at(i, pos);
	}
	//ノットベクトルを求めてオフセット曲線を生成する
	int error;
	ofs_brep = MGLBRep(dataPoint, bp1, knotVector, error);	//精度十分の曲線を生成する
	if(error)return false;
	return true;
}

//法線ベクトルが指定されていないC2連続曲線のオフセット
int MGCurve::offset_c2_proc(
	const MGLBRep& ofs_value_lb,	//オフセット量
	MGLBRep& ofs_brep,			//オフセットカーブ
	MGUnit_vector& preN,		//前回のノーマルベクトル
	int& freverse			//向きを逆にしているというフラグ
)const{
	MGKnotVector knotVector = offset_make_knotvector(ofs_value_lb);	//十分分割したノットベクトルを求める
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);

	//制御点を生成する
	int len = dataPoint.length();
	MGBPointSeq bp1(len, sdim());
	for(int i = 0; i < len; i++){
		MGUnit_vector T, N, B;
		double curvature, torsion, ofs_value = ofs_value_lb.eval_position(dataPoint(i)).ref(0);
		Frenet_frame(dataPoint(i), T, N, B, curvature, torsion);
		if(MGAZero(curvature))N = preN;		//曲率が小さいときノーマルは前のを使う
		//ノーマルが全部同じ方向になるようにする(180度以上開かない)
		if((preN % N) < 0)freverse = !(freverse);
		int fdir = -2 * freverse + 1;	//freverseがfalseのとき1, trueのとき-1
		if(!MGMZero(curvature))	//ゼロ割チェック
			if(ofs_value * fdir > (1. / curvature))return false;	//曲率半径より大きいオフセットは認めない
		preN = N;	//ノーマルを取っておく
		MGPosition pos;
		pos = eval(dataPoint(i));
		pos += fdir * N * ofs_value;
		bp1.store_at(i, pos);
	}

	//ノットベクトルを求めてオフセット曲線を生成する
	int error;
	ofs_brep = MGLBRep(dataPoint, bp1, knotVector, error);	//精度十分の曲線を生成する
	if(error)return false;
	return true;
}

//曲線を折れで分割する(オフセット量を示す曲線の折れついても分割する)
int MGCurve::divide_multi_ofs(
	const MGLBRep& ofs_value_lb,		//オフセット量を示す曲線
	MGPvector<MGCurve>& crv_list		//分割した曲線リスト
)const{
	crv_list.clear();
	//オフセット量を示す曲線と元の曲線のパラメータレンジが違うとエラー
	if(param_range() != ofs_value_lb.param_range())return 0;
	size_t	start_index = ofs_value_lb.order() - 1, index = 0, count = 0, multi = 0, vbdim = ofs_value_lb.bdim();
	MGKnotVector value_knot = ofs_value_lb.knot_vector();
	do{
		if(ofs_value_lb.order() == 2){	//オーダー2のときの処理
			index = start_index + 1; multi = 1;
		}else{
			multi = value_knot.locate_multi(start_index, 2, index);
		}
		MGCurve *temp_crv = part(value_knot(start_index), value_knot(index));
		count += temp_crv->divide_multi(crv_list);
		delete temp_crv;
		start_index = index + multi - 1;
	}while(index != vbdim && value_knot(start_index) < param_e());	//多重度が見つからなかったら終わり
	return count;
}

//1スパンの分割数を求める
int MGCurve::offset_div_num(
	const MGInterval& interval	//分割数を求めるパラメータ範囲
)const{
	//２＊order()で分割して最大2次微分値を求める
	int ord = order();
	if(ord<=2) ord = 4;		//Ellipse, Straightのとき
	int ord2=ord*2;
	double max_deriv = 0.0, tpara = interval.low_point();
	double shortspan = (interval.high_point() - tpara)/ord2;
	double oneSpanLength = interval.length().value();
	double oneSpanLength2=oneSpanLength*oneSpanLength;
	for(int j = 0; j <=ord2; j++){
		//パラメータ範囲(0,1)の2次微分値に直すためにパラメータ範囲の2乗をかける
		double deriv = eval(tpara, 2).len() * oneSpanLength2;
		if(deriv > max_deriv) max_deriv = deriv;
		tpara += shortspan;
	}

	double onelzero=1./MGTolerance::line_zero();
	//分割数は (1 / tol)^POW * sqrt(max_deriv / 8)で最小値は order()である
	double sqrt_max_deriv8=sqrt(max_deriv / 8.);
	int ndiv = int(pow(onelzero, POW_HIGH) * sqrt_max_deriv8);
	if(ndiv > NUM_DIV)
		ndiv = int(pow(onelzero, POW_LOW) * sqrt_max_deriv8);
	if(ndiv < ord2) ndiv = ord2;
	else if(ndiv > NUM_DIV*2) ndiv=NUM_DIV*2;
	return ndiv;
}

//向きが同じB表現曲線リストを接続する(LBRepのみ)。join_crvlに接続した曲線リストが入る。
//戻り値は、引数の曲線リストの向きが違うとき、同じB表現同士でなかったときfalseが返る。
int join(
	MGPvector<MGCurve>& crvl,
	MGPvector<MGCurve>& join_crvl
){
	int num = crvl.size(), rc = 0;
	MGCurve *pre_pcrv = 0, *cur_pcrv, *next_pcrv;
	if(!num)
		return false;

	if(num == 1){	//曲線が１本のときの処理
		cur_pcrv = crvl[0]->clone();
		cur_pcrv->remove_knot();
		join_crvl.push_back(cur_pcrv);
		return 1;
	}

	cur_pcrv = crvl[0]->clone();
	for(int i = 0; i < num - 1; i++){
		next_pcrv = crvl[i+1];
		if(!cur_pcrv || !next_pcrv)
			return false;

		MGLBRep* lb1=dynamic_cast<MGLBRep*>(cur_pcrv);
		MGLBRep* lb2=dynamic_cast<MGLBRep*>(next_pcrv);
		if(!lb1 || !lb2)
			return false;

		pre_pcrv = join2LBRep(*lb1,*lb2);
		if(!pre_pcrv){
			rc++;
			cur_pcrv->remove_knot();
			join_crvl.push_back(cur_pcrv);
			delete pre_pcrv;
			cur_pcrv = next_pcrv->clone();
			continue;
		}
		delete cur_pcrv;
		cur_pcrv = pre_pcrv;
	}
	cur_pcrv->remove_knot();
	join_crvl.push_back(cur_pcrv);
	rc++;
	return rc;
}

/*//向きが同じ2本のB表現曲線を接続する(同じ種類のとき)
MGCurve* MGCurve::join(const MGCurve& crv1)const{
	return NULL;
}*/

/*//向きが同じ2本のB表現曲線を接続する(同じ種類のとき)
MGCurve* MGRLBRep::join(const MGCurve& crv1) const{
	const MGRLBRep *prlbrep1 = dynamic_cast<const MGRLBRep*>(&crv1);
	if(!prlbrep1)return NULL;
	int which, cont;
	double ratio;
	cont = this->continuity(*prlbrep1, which, ratio);
	if(cont < 0 || which == 0 || which == 3)return NULL;
	return new MGRLBRep(*this, cont, which, *prlbrep1);
}*/

//曲線をオフセットするのに十分分割したノットベクトルを返却する
//オフセット量曲線も考慮に入れ、分割数の多い方にあわせている。
MGKnotVector MGCurve::offset_make_knotvector(
	const MGLBRep& ofs_value_lb
)const{
	//元となるノットベクトルを生成する
	const MGKnotVector& tempKnotVector = knot_vector();
	MGKnotVector knotVector(tempKnotVector, 4);	//オーダー4のノットベクトルに作り替える
	for(size_t i=tempKnotVector.order()-1; i<tempKnotVector.bdim(); i++){
		double	tpara = 0.0,					//テンポラリ
				spara = tempKnotVector(i),		//スパンの始点
				epara = tempKnotVector(i + 1);	//スパンの終点
		if(epara - spara < param_error())
			continue;	//マルチノットのときの処理(RLBRepのみ)

		//1スパンの分割数を決定する
		MGInterval interval(spara, epara);
		size_t ndiv = offset_div_num(interval),					//オフセット曲線の分割数
			tmp_ndiv = ofs_value_lb.offset_div_num(interval);	//オフセット量曲線の分割数
		if(tmp_ndiv>ndiv)
			ndiv=tmp_ndiv;						//分割数の多い方を用いる
		double shortspan = (epara - spara) / ndiv;
		tpara = spara;
		for(size_t j=0; j<ndiv; j++){
			knotVector.add_data(tpara);
			tpara += shortspan;
		}
	}
	return knotVector;
}
