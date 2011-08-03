/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Bisection.h"
#include "mg/Box.h"
#include "mg/Knot.h"
#include "mg/Unit_vector.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/SurfCurve.h"
#include "mg/FSurface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;
#define EXTEND_COEF 2.0		//ベクトル投影のスイープ面の長さを求めるのに使用

// Implementation of MGFSurface projection.

//向きが同じ2本のB表現曲線を接続する(同じ種類のとき)
MGLBRep* join2LBRep(const MGLBRep& crv1, const MGLBRep& crv2);

//向きが同じ５次(xyzuv)のB表現曲線リストを接続し、xyzを見てline_zeroでリムーブノットする
//join_crvl_uv,join_crv_xyzに接続した曲線リストが入る。
//Function's return value is the number of output curves.
void prjJoin(
	MGPvector<MGLBRep>& crvl,
	MGPvector<MGCurve>& join_crvl_uv,
	MGPvector<MGCurve>& join_crvl_xyz
){
	int num = crvl.size();
	if(!num)
		return;

	if(num==1){	//曲線が１本のときの処理
		const MGLBRep& lbrep=*(crvl[0]);
		join_crvl_uv.push_back(new MGLBRep(2,lbrep,0,3));
		join_crvl_xyz.push_back(new MGLBRep(3,lbrep,0,0));
		return;
	}

	MGLBRep *cur_pcrv;
	MGLBRep *next_pcrv;
	cur_pcrv = crvl.release(0);
	for(int i=1; i<num; i++){
		next_pcrv = crvl.release(i);
		MGLBRep *pre_pcrv = join2LBRep(*cur_pcrv,*next_pcrv);
		if(pre_pcrv){//If successfully joined.
			delete cur_pcrv;
			delete next_pcrv;
			cur_pcrv = pre_pcrv;
		}else{
			join_crvl_uv.push_back(new MGLBRep(2,(*cur_pcrv),0,3));
			join_crvl_xyz.push_back(new MGLBRep(3,(*cur_pcrv),0,0));
			delete cur_pcrv;
			cur_pcrv = next_pcrv;
		}
	}
	join_crvl_uv.push_back(new MGLBRep(2,*cur_pcrv,0,3));
	join_crvl_xyz.push_back(new MGLBRep(3,*cur_pcrv,0,0));
	delete cur_pcrv;
}

//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int MGFSurface::project(
	const MGCurve& crv,
	MGPvector<MGCurve>& vec_crv_uv,	
		//Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	const MGCompositeCurve* ccrv=dynamic_cast<const MGCompositeCurve*>(&crv);
	if(ccrv){
		ccrv->project_onto_surface(*this,vec_crv_uv,vec_crv,vec);
		return vec_crv.size();
	}
	//初期化
	assert(crv.sdim() > 1);
	std::auto_ptr<MGCurve> crvRemoveKnot(crv.clone());
	crvRemoveKnot->remove_knot();

	//ベクトル投影かどうか調べる
	if(!vec.is_null())	//ベクトル投影
		return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, vec);

	const MGPlane* pln=dynamic_cast<const MGPlane*>(this);
	if(pln)
		//Planeでベクトルが指定されていないときはplaneのノーマルでベクトル投影する
		return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, pln->normal());

	MGPvector<MGCurve> divided_curves;
	//マルチノットになるときは切って投影する
	int n=crvRemoveKnot->divide_multi(divided_curves);
	MGPvector<MGLBRep> crv_xyzuv_vector;//Projected curve of(x,y,z,u,v).
	for(int i=0; i<n; i++){
		project1span(*(divided_curves[i]),crv_xyzuv_vector);
	}

	//投影曲線列の冗長ノットを削除して接続する
	prjJoin(crv_xyzuv_vector, vec_crv_uv, vec_crv);
	return vec_crv.size();
}

//Obtain the projected curve of a curve onto the surface.
//project1span does not divide the curve at the c0 continuity, which is treated by
//project(). See projec().
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project1span' is curves of space dimension 5, which are (x,y,z,u,v).
//Here (x,y,z) is the world coordinates, and (u,v) is parameter of this surface.
void MGFSurface::project1span(
	const MGCurve& crv,	//The target curve to project.
	MGPvector<MGLBRep>& crv_xyzuv_vector//Projected curve of(x,y,z,u,v) will be appended.						
)const{
	//パラメータ範囲を求める
	std::deque<MGPosition> ranges;//ranges[2*i] to ranges[2*i+1] is the projection range.
			//Let ranges[.]=tuv, then tuv[0] is the parameter of the curve, and (tuv[1], tuv[2])
			//is the parameter of this surface (u,v).

	double sparam = crv.param_s();
	int counter=0, rc;
	do{
		MGPosition tuv[2];
		rc = prj2GetParamRange(crv,counter,tuv,counter);
		if(rc==1 || rc==-1){
			ranges.push_back(tuv[0]);
			ranges.push_back(tuv[1]);
		}
	}while(rc < 0);	//crv1の途中であれば繰り返す

	//投影を行う
	while(ranges.size()>=2){
		//ranges will be pop_front by prj2OneCurve.
		MGLBRep* projectedCurve;
		prj2OneCurve(crv,ranges,projectedCurve);

		//投影曲線列の冗長ノットを削除して接続する
		if(projectedCurve)
			crv_xyzuv_vector.push_back(projectedCurve);
	}
}

//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is general world coordinate curves('vec_crv')
int MGFSurface::project(
	const MGCurve& crv,			//given curve.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec			//projection vector.
		//if vec = NULL then calculate perpendicular project.
)const{
	MGPvector<MGCurve> vec_crv_uv;
	return project(crv,vec_crv_uv,vec_crv,vec);
}

//面直に投影した点を返却する
//戻り値は、交点または面直点が求まったときは1、求まらなかったときは0を返却する
int MGFSurface::project_normal(
	const MGPosition& pos,
	const MGPosition& uv_guess,	//推量パラメータ
	MGPosition& uv
)const{
	double dPreTol = MGTolerance::set_wc_zero(MGTolerance::wc_zero()*0.5);
	int rc = perp_point(pos, uv, &uv_guess);
	MGTolerance::set_wc_zero(dPreTol);
	return rc;
}

//ベクトル投影は、カーブを折れで分割して行い、後で接続する
int MGFSurface::projVector(
	const MGCurve& crv,
	MGPvector<MGCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	MGPvector<MGCurve> mult_crvl;
	MGPvector<MGLBRep> temp_vec_crv;
	//マルチノットになるときは切って投影して再度接続する
	int num_mult = crv.divide_multi(mult_crvl);
	for(int i=0; i<num_mult; i++)
		projVectorProc(*(mult_crvl[i]), temp_vec_crv, vec);
	prjJoin(temp_vec_crv, vec_crv_uv, vec_crv);
	return vec_crv.size();
}

//ベクトル投影は、crvをvec方向にスイープした面と元の面との交線で求める
int MGFSurface::projVectorProc(
	const MGCurve& crv,
	MGPvector<MGLBRep>& crv_xyzuv_vector,//Projected curve of(x,y,z,u,v) will be appended.						
	const MGVector& vec
)const{
	MGBox box(this->get_box() | crv.box());
	double extLeng = box.len() * EXTEND_COEF;	//安全の為
	const double srfError=param_error();	
	//おそらくplaneの時の処理を別に書いてやらないとうまく行かないと思う
	//精度的に問題なければこの処理をコメントにしても大丈夫かも2001/11/26
	if(get_surface_pointer()->type() == MGSURFACE_PLANE){
		//曲線の中心からベクトル方向の直線を生成する
		MGStraight centerStraight(MGSTRAIGHT_UNLIMIT,vec,crv.center());
		//直線と面との交点を求める
		MGCSisect_list isectList = isect(centerStraight);
		if(isectList.size() > 0){
			std::list<MGCSisect>::const_iterator iterIsect = isectList.begin();
			MGPosition posIsect = iterIsect->point();
			//延長する長さを直線距離+曲線長さにする
			extLeng = (posIsect - crv.center()).len() + crv.length();
			extLeng *= EXTEND_COEF;	//安全の為
		}
	}
	if(MGAZero(extLeng / EXTEND_COEF))extLeng = 1.0;		//絶対ゼロ以下の場合の処理
	std::auto_ptr<MGSurface> apCutSrf(crv.sweep(vec, extLeng, -extLeng));
	if(!apCutSrf.get())
		return 0;
	double errsave=MGTolerance::set_wc_zero(MGTolerance::line_zero());
	MGSSisect_list isectl = isect(*apCutSrf);//cout<<isectl<<endl;/////////////********
	MGTolerance::set_wc_zero(errsave);
	int i;
	MGSSisect_list::SSiterator ssiter;
	for(ssiter = isectl.begin(), i = 0; ssiter != isectl.end(); ssiter++, i++){
		//カーブの向きをあわせてからMGPvectorにpushする
		std::auto_ptr<MGCurve> apCrv2d(ssiter->release_param1());
		std::auto_ptr<MGCurve> apCrv3d(ssiter->release_line());
		MGCurve &crvOth2d = ssiter->param2();
		assert(apCrv2d.get() && apCrv3d.get());
		//uパラメータが小さい方を始点とする(元カーブの向きと同じにする)
		double spt = (crvOth2d.start_point())(0), ept = (crvOth2d.end_point())(0);
		if(spt > ept){
			apCrv2d->negate();
			apCrv3d->negate();
		}
		//cout<<*apCrv2d<<*apCrv3d<<endl;
		double	param_len = apCrv2d->box().len(),	//ショートアークのチェックをする
				world_len = apCrv3d->box().len();
		if(param_len < srfError || world_len < MGTolerance::wc_zero())continue;
		const MGLBRep *pCrv2d = dynamic_cast<const MGLBRep*>(apCrv2d.get());
		const MGLBRep *pCrv3d = dynamic_cast<const MGLBRep*>(apCrv3d.get());
		assert(pCrv2d && pCrv3d);
		size_t bdim = pCrv2d->bdim();
		MGBPointSeq bpXYZUV(5,pCrv3d->line_bcoef());
		const MGBPointSeq &bpUV = pCrv2d->line_bcoef();
		for(size_t i2=0; i2<bdim; i2++)
			bpXYZUV.store_at(i2,bpUV(i2),3,0);
		//cout<<bpXYZUV<<endl;
		crv_xyzuv_vector.push_back(new MGLBRep(pCrv2d->knot_vector(), bpXYZUV));
	}
	return i;
}

//トレランスに応じた細分したデータポイントを求める
void MGFSurface::prjGetDataPoints(
	const MGCurve& curve,
	const MGPosition& tuv0,//parameter range of curve to get the data point,
	const MGPosition& tuv1,//From tuv0 to tuv1.
	MGNDDArray& tau		//data points will be output.
)const{
	double t0=tuv0[0], t1=tuv1[0];
	MGInterval cspan(t0, t1);
	int n1=curve.offset_div_num(cspan);

	const MGSurface* srf=get_surface_pointer();

	MGCurve* pcrv1=srf->parameter_curve(1,(tuv0[1]+tuv1[1])*.5);
	MGInterval cspan2(tuv0[1], tuv1[1]);
	int n2=pcrv1->offset_div_num(cspan2);
	delete pcrv1;

	MGCurve* pcrv2=srf->parameter_curve(0,(tuv0[2]+tuv1[2])*.5);
	MGInterval cspan3(tuv0[2], tuv1[2]);
	int n3=pcrv2->offset_div_num(cspan3);
	delete pcrv2;

	int n=n1;
	if(n<n2)
		n=n2;
	if(n<n3)
		n=n3;

	tau.resize(n+1);
	double delta=(t1-t0)/n;
	double t=t0;
	for(int i=0; i<n; i++, t+=delta){
		tau(i)=t;
	}
	tau(n)=t1;
}

//Get world position data at (u, v) and store the coordinates
//and the parameter in bp.
void MGFSurface::prj_store_bp(
	int i,
	double u, double v,
	MGBPointSeq& bp
)const{
	MGPosition Psrf=eval(u, v);	//始点を求める
	bp.store_at(i,Psrf,0,0,3);
	bp(i,3)=u; bp(i,4)=v;
}

//投影を実行して５次元曲線を取得する
//curveにおけるデータポイントと投影ベクトルから投影点列を作成し
//点列から５次元(xyz,uv)の曲線列を作成して返却する
void MGFSurface::prj2OneCurve(
	const MGCurve& curve,	//target curve to prject.
	std::deque<MGPosition>& ranges,//start(ranges[0]) and end(ranges[1]) point parameter
							//of the curve and the face.
							//On return ranges will be so updated that processed rages[0] to [1]
							//are pop_front().
	MGLBRep*& crvProjected	//newed object will be returend when obtained.
							//When not obtained, null will be returned.
)const{
	crvProjected=0;
	MGPosition tuv0=ranges.front(); ranges.pop_front();
	MGPosition tuv1=ranges.front(); ranges.pop_front();

	MGNDDArray tau;//data point to get the projected points.
				//tau[0] is tuv[0][0], and tau[n-1]=tuv[1][0].
	prjGetDataPoints(curve,tuv0,tuv1,tau);

	int ntau = tau.length();
	int ntaum1=ntau-1;

	MGBPointSeq bp_add(ntau, 5);	//xyz,uvをあわせたＢ表現[x,y,z,u,v]
	prj_store_bp(0,tuv0[1],tuv0[2],bp_add);//始点を求める

	const MGSurface* srf=get_surface_pointer();
	MGPosition uv, uv_guess(tuv0[1],tuv0[2]);
	int i=1;
	for(; i<ntaum1 ; i++){
		MGPosition Pcrv = curve.eval_position(tau(i));
		project_normal(Pcrv, uv_guess, uv);
		//	MGPosition tuv2=proj2GetParamIter(curve,tuv1,tau[i],tuv1[0]);
		//	ranges.push_front(tuv1);//to process in the next call of prjOneCurve.
		//	ranges.push_front(tuv2);//Next process's start .
		//	MGPosition tuv(tau[i-1],uv_guess[0],uv_guess[1]);
		//	tuv1=proj2GetParamIter(curve,tuv,tau[i-1],tau[i]);//New end.
		//	break;
		
		size_t periNum;
		srf->on_a_perimeter(uv(0), uv(1), periNum);
			//辺にのっていたらそのパラメータを使用する(uvがupdateされている）
		uv_guess = uv;		//推量点を更新する
		prj_store_bp(i,uv[0],uv[1],bp_add);//始点を求める
	}
	prj_store_bp(i++,tuv1[1],tuv1[2],bp_add);//End pointを求める
	bp_add.set_length(i);

	//曲線を作成してベクトルに挿入する
	MGBox bxAddXYZ(3,bp_add.box());	//xyzのボックスを求める
	//ショートアークのチェックをする
	if(bxAddXYZ.len() <= MGTolerance::wc_zero())
		return;

	MGBPointSeq tmpBpXYZ(3,bp_add);
	MGNDDArray tauXYZ(tmpBpXYZ);
	crvProjected = new MGLBRep(tauXYZ, bp_add, 4, 4.);
	crvProjected->remove_knot(0,3);	//最初の３次元(XYZ)をline_zeroでノット削除
		//cout<<(*crvProjected);
}

//Normalize the parameter of the input range[.][0], i.e. if the value is alomost equalt to
//a knot value, round the value into the knot value. This is usefull especially for
//start or end parameter value.
int normalize_check_isect(
	const MGCurve& curve,
	int retval,
	MGPosition range[2]
){
	//パラメータ値を正規化する
	MGPosition& tuv0=range[0];
	MGPosition& tuv1=range[1];
	double ts=tuv0(0)=curve.param_normalize(tuv0[0]);
	double te=tuv1(0)=curve.param_normalize(tuv1[0]);
	if((te-ts)<=curve.param_error()){
		if(retval==1)
			return 0;	//投影範囲が誤差範囲の場合交点となるので0を返す
		else
			return -2;
	}
	return retval;
}

//スタートパラメータを与え投影可能なパラメータ範囲を1つだけ取得する。
//戻り値は	1:範囲が求まった(thisの最後まで)
//			0:範囲が求まらなかった(thisの最後まで)
//			-1:範囲が求まった(thisの途中まで、まだ面直範囲があるかもしれない)
//			-2;範囲が求まらなかったが、thisの途中で、まだ面直範囲があるかもしれない
int MGFSurface::prj2GetParamRange(
	const MGCurve& curve,
	int start_counter,	//input start counter of the curve parameter incrementation.
	MGPosition range[2],
		//range[0][0]=start point parameter of curve,
		//range[0][1-2]=the corresponding surface(this surface's) parameter (u,v)
		//of the start point of the range.
		//Regarding to range[1] : the same.
	int& next_counter	//Updated curve parameter incremental counter will be output.
)const{
	//初期化
	//std::cout<<curve<<std::endl;
	int ndiv = curve.intersect_dnum()+2;

	//curveを分割し、投影可能か調べる
	const double delta=curve.param_span()/ndiv;	//チェックポイントのパラメータスパン
	const double tend=curve.param_e();

	int t_is_on_surface;
	double t=curve.param_s()+delta*double(start_counter), t_pre;
	MGPosition uv;
	int i=start_counter;//counter for t's incremental number.
	if(t_is_on_surface=perp_one(curve.eval(t),uv)){
		range[0]=MGPosition(t,uv[0],uv[1]);
	}else{
		//Find the 1st t that is on surface.
		while(!t_is_on_surface && i<ndiv){
			t_pre=t;
			t+=delta; i++;
			t_is_on_surface=perp_one(curve.eval(t),uv);
		}
		if(t_is_on_surface){
			range[0]=proj2GetParamIter(curve,MGPosition(t,uv[0],uv[1]),t_pre,t);
		}else{
			if(perp_one(curve.eval(tend),uv)){
				range[0]=proj2GetParamIter(curve,MGPosition(tend,uv[0],uv[1]),tend-delta,tend);
				range[1]=MGPosition(tend,uv[0],uv[1]);
				return normalize_check_isect(curve,1,range);
			}
			return 0;
		}
	}

	//Here t_is_on_surface=true always holds.
	MGPosition uv_save;
	while(t_is_on_surface && t<tend && i<ndiv){
		uv_save=uv;
		t_pre=t;
		t+=delta; i++;
		t_is_on_surface=perp_point(curve.eval(t),uv,&uv_save);
	}

	int retval;
	if(!t_is_on_surface){
		range[1]=proj2GetParamIter(curve,MGPosition(t_pre,uv_save[0],uv_save[1]),t_pre,t);
		retval=-1;
	}else{
		uv_save=uv;
		if(perp_point(curve.eval(tend),uv,&uv_save)){
			range[1]=MGPosition(tend,uv[0],uv[1]);
		}else{
			range[1]=proj2GetParamIter(curve,MGPosition(t,uv[0],uv[1]),t,tend);
		}
		retval=1;
	}

	//終了処理
	next_counter=i;		//次回の検索開始パラメータを返す
	return normalize_check_isect(curve,retval,range);
}

class MGPrjBisect{
public:

MGPrjBisect(
	double ts,	//parameter range from ts to te.
	double te
);

//Virtual Destructor
virtual ~MGPrjBisect(){;};

//compare with the previous function value(the initial value is set
//by set_initial_t) and replace t with the previous one if necessary.
//The function's return value is the new parameter value.
virtual void compare_replace(
	double t	//parameter value to compare at.
)=0;

int solve(
	double tolerance//The tolerance to halt the bisection iteration.
);

protected:
	double m_ts, m_te;
		//the curve's parameter range from m_ts to m_te, or from m_te to m_ts.
		//note that m_ts<m_te does not always hold.
};

//A virtual super class to solve non-linear equations by bicection methos.
MGPrjBisect::MGPrjBisect(
	double ts,	//parameter range from ts to te.
	double te
):m_ts(ts), m_te(te){}

//Compute the fn(t)'s parameter value that is the maxima 
//Function's return value will be the solution obtained.
int MGPrjBisect::solve(
	double tolerance//The tolerance to halt the bisection iteration.
){
	assert(tolerance>0.);
	int nrepition=0;
	double span=m_te-m_ts;
	tolerance+=MGTolerance::mach_zero()*2.;
	while(fabs(span)>=tolerance){//m_ts > m_te may happend.
		span*=.5; nrepition++;
	    double t=m_ts+span;
		compare_replace(t);
	}
	return nrepition;
}

class MGProjBoundaryParam: public MGPrjBisect{
public:
	MGProjBoundaryParam(
		const MGFSurface& surf,
		const MGCurve& curve,
		const MGPosition& tuv,
		double ts, double te
	):MGPrjBisect(ts,te),m_surf(surf),m_curve(curve),m_uv(tuv[1],tuv[2]){
		if(tuv[0]==ts){
			m_te=ts;
			m_ts=te;
		}
	//m_te holds the curve parameter value whose point has
	//a perpendicular point ont the surface.
	}

	//compare with the previous function value(the initial value is set
	//by set_initial_t) and replace t with the previous one if necessary.
	//The function's return value is the new parameter value.
	virtual void compare_replace(
		double t	//parameter value to compare at.
	){
		MGPosition uv;
		MGPosition P=m_curve.eval(t);
		if(m_surf.project_normal(P,m_uv, uv)){
			m_uv=uv;
			m_te=t;
		}else{
			m_ts=t;
		}
	}

	MGPosition get_result(){return MGPosition(m_te,m_uv[0],m_uv[1]);};
	const MGFSurface& m_surf;
	const MGCurve& m_curve;
	MGPosition m_uv;
};

//curve上の投影可能なパラメータと投影不可能なパラメータを与えて、
//間にある投影可能な境界パラメータ値を求めて返却する
//チェックポイントの移動量 < (パラメータ範囲*rc_zero())になれば終了。
//Function's return value is MGPosition tuv, where
//tuv[0]=start point parameter of curve,
//(tuv[1], tuv[2])=the corresponding surface(this surface's) parameter (u,v)
//of the start point of the range.
MGPosition MGFSurface::proj2GetParamIter(
	const MGCurve& curve,//target curve to project.
	const MGPosition& tuv,//tuv[0] is the curve's parameter value that has a perp point onto the
				//surface. and (tuv[1], tuv[2]) is the surface's parameter of the perp point.
				//tuv[0] is eithere ts or te.
	double ts,	//(ts, te) is curve's parameter range that indicates
	double te	//between ts and te there must be the boundary point.
				//That is one of the following situations occurs:
				//(1) there are no perp points at ts and there is a perp point at te,
				//(2) there is a perp point at ts and there are no perp points at te,
)const{
	//初期化
	const double crvError = curve.param_error()*2.;//, srfError = get_surface_pointer()->param_error();
	MGProjBoundaryParam solver(*this,curve,tuv,ts,te);
	solver.solve(crvError);
	return solver.get_result();
}
