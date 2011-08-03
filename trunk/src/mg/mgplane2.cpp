/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Interval.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Curve.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/SSisect_list.h"
#include "mg/Point.h"
#include "mg/Position_list.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Bvi2pl.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGPlane.cc
//
// MGPlaneクラスは３次元空間における平面を表すクラスである。
// MGPlaneクラスでは以下のようなパラメータ表現を使用する。
// Point = m_root_point + u * m_uderiv + v * m_vderiv
using namespace std;

////// Member function(intersection related ones)メンバ関数 //////

// Surface と Curve の交点を求める。
//Compute intesection of Plane and Curve.
MGCSisect_list MGPlane::isect(const MGCurve& curve) const{
	return curve.isect(*this);
}

MGCSisect_list MGPlane::isect(const MGStraight& curve) const{
	MGCSisect_list list(&curve,this);

	const MGVector& lvec=curve.direction();
	const MGPosition& lpoint=curve.root_point();
	double r=lvec%m_normal;
	if(!MGMZero(r)){
		double t=(m_d-lpoint%m_normal)/r;
		if(curve.in_range(t)){
			MGPosition ipoint=lpoint+lvec*t;
			list.append(MGCSisect(ipoint,t,uv(ipoint)));
		}
	}
	return list;
}

MGCSisect_list MGPlane::isect(const MGRLBRep& curve) const{
	return curve.isect(*this);
}

MGCSisect_list MGPlane::isect(const MGEllipse& el) const{
	MGCSisect_list list(&el,this);
	if(!has_common(el)) return list;

	if(!m_normal.parallel(el.normal())){
		MGStraight line; MGCCisect is;
		//Compute intersection line with the plane that includes el.
		relation(MGPlane(el.normal(),el.center()), line);
			//line is the line of the intersection of two planes.
		MGCCisect_list list2=el.isect(line);
			//list2 is the intersection points list of line and ellipse.
		size_t n=list2.entries();
		for(size_t i=0; i<n; i++){
			is=list2.removeFirst();
			MGPosition ipoint(is.point());
			list.append(MGCSisect(ipoint,is.param1(),uv(ipoint)));
		}
	}
	return list;
}

MGCSisect_list MGPlane::isect(const MGLBRep& curve) const{
	return curve.isect(*this);
}

MGCSisect_list MGPlane::isect(const MGSurfCurve& curve) const{
	return curve.isect(*this);
}

MGCSisect_list MGPlane::isect(const MGBSumCurve& curve) const{
	return curve.isect(*this);
}

//Compute intersections with MGLBRep lb that does not have C0 continuity in it.
MGCSisect_list MGPlane::isect_withC1LB(const MGLBRep& lb)const{
	return lb.intersect_with_plane(*this);
}

//isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisect_list MGPlane::isect_with_noCompoSC(const MGSurfCurve& scrv)const{
	return scrv.isect_noCompo(*this);
}

// Plane と Surface の交線を求める。
// Compute intesection of Plane and Surface.
MGSSisect_list MGPlane::isect(const MGSurface& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGPlane::isect(const MGPlane& srf2) const{
	MGSSisect_list list(this, &srf2);
	double g1[4], g2[4], line[2][3]; int iflag;
	for(size_t i=0; i<3; i++){
		g1[i]=m_normal.ref(i); g2[i]=srf2.m_normal.ref(i);
	}
	g1[3]=m_d; g2[3]=srf2.m_d;
	bvi2pl_(g1,g2,line[0],&iflag);	//Compute intersection line.
	if(iflag==1){
		MGVector dir(3,line[1]); MGPosition p(3,line[0]);
		MGStraight* sl=new MGStraight(MGSTRAIGHT_UNLIMIT, dir, p);
		MGStraight* param1=new MGStraight;
		param1->set_straight(MGSTRAIGHT_UNLIMIT,uv(dir),uv(p));
		MGStraight* param2=new MGStraight;
		param2->set_straight(MGSTRAIGHT_UNLIMIT,srf2.uv(dir),srf2.uv(p));
		MGSSisect ssi(sl,param1,param2,MGSSREL_ISECT);
		list.append(ssi);
	}
	return list;
}
MGSSisect_list MGPlane::isect(const MGSphere& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGPlane::isect(const MGCylinder& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGPlane::isect(const MGSBRep& srf2) const{
	MGSSisect_list list=srf2.intersectPl(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGPlane::isect(const MGRSBRep& srf2) const{
	MGSSisect_list list=srf2.intersectPl(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGPlane::isect(const MGBSumSurf& srf2) const{
	MGSSisect_list list=srf2.intersectPl(*this);
	list.replace12();
	return list;
}

// "isect_guess_straight" computes one intersection point of surface and
//a straight line, given initail guess parameter values of the surface and 
//the straight line.
//Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int MGPlane::isect_guess_straight(
	const MGStraight& sl,	//Straight line.
	double ti,			//Initial guess parameter value of the straight.
	const MGPosition& uvi,	//Input initial guess parameter value
						// of the i.p. of the surface. 
	double& t,			//Straight parameter obtained.
	MGPosition& uvout	//Surface parameter value obtained(u,v). 
)const{
	const MGVector& lvec=sl.direction();
	const MGPosition& lpoint=sl.root_point();
	double r=lvec%m_normal;
	if(MGMZero(r)) return 0;
	t=(m_d-lpoint%m_normal)/r;
	MGPosition ipoint=lpoint+lvec*t;
	uvout=uv(ipoint);
	return sl.in_range(t);
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr
//(See MGSurface), will get shortest parameter line
//necessary to compute intersection.
MGCurve* MGPlane::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	size_t incr			//Incremental valuse of B-coef's id.
) const{
	u=uv.ref(0); v=uv.ref(1);
	MGVector direction;
	MGPosition P;
	if(kdt){ u+=du; direction=v_deriv(); P=eval(u,0.);}
	else   { v+=dv; direction=u_deriv(); P=eval(0.,v);}
	MGStraight* sl=new MGStraight();
	sl->set_straight(MGSTRAIGHT_UNLIMIT,direction,P);
	return sl;
}

//"isect_inner_dt" is a dedicated function of isect_start,
// comutes adequate incremental parameter value and parameter line kind
//(u=const or v=const).
void MGPlane::isect_inner_dt(
	size_t n,				//num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,//intersection point obtained last(of this).
	double& du, double& dv,	//incremental length from previous to uvnow is input.
				//New du or dv will be output according to kdt's return value.
	int& kdt,	//Parameter kind used so far is input, will be output as:
				//=1:parameter line kind(u=const), =0: v=const,
				//=-1:should halt computation since incremental value is zero.
	double acuRatio	//Accurate ratio.
)const{
	const MGVector& dfdu=u_deriv();
	const MGVector& dfdv=v_deriv();
	double ulen=dfdu.len(), vlen=dfdv.len();
	double error=MGTolerance::rc_zero();
	double uerr=ulen*error, verr=vlen*error;
	double abdu=fabs(du), abdv=fabs(dv);
	double ratio=acuRatio;
	if(abdu<uerr){
		if(abdv<verr){kdt=-1; return;}
		else kdt=0;
	}else if(abdv<verr) kdt=1;
	else{
		double fuu=ulen*ulen, fuv=dfdu%dfdv, fvv=vlen*vlen;
		double dffu=fuu*du+fuv*dv, dffv=fuv*du+fvv*dv;
		double dffubyfvv=dffu*dffu*fvv, dffvbyfuu=dffv*dffv*fuu;
		if(kdt==-1){
			if(dffubyfvv>=dffvbyfuu) kdt=1; else kdt=0;
		}else if(dffubyfvv>=dffvbyfuu*2.){
			//u=const and v-varying parameter line.
			if(abdu>uerr*4.) kdt=1;
		}else if(dffvbyfuu>=dffubyfvv*2.){
			//v=const and u-varying parameter line.
			if(abdv>verr*4.) kdt=0;
		}
		if(kdt) ratio*=dffubyfvv; else ratio*=dffvbyfuu;
		ratio/=(dffubyfvv+dffvbyfuu);
	}
//Define new dt,kdt.
	double dt=4., dtold;
	dt*=isect_dt_coef(n);
	if(kdt){ dt*=ulen; dtold=du;}
	else   { dt*=vlen; dtold=dv;}
	if(dtold<0.) dt=-dt;
	dt*=isect_dt_coef(n)*ratio;
	if(n){
		//When this is not the 1st call of isect_inner_dt,
		//dt must not exceed twice or half of the old dt.
		double dtr1=dt*dt, dtr2=dtold*dtold;
		if(dtr1 > 2.*dtr2) dt=dtold*2.;
		else if(dtr1 < .5*dtr2) dt=dtold*.5;
	}
	if(kdt) du=dt; else dv=dt;
	return;
}

//isect_startHPL compute one intersection line of two surfaces, this and sf2,
// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2).
// isect_startHPL halts the computation when intersection
// reached to a boundary of this or sf2, or reached to one of the points
// in uvuv_list.
//The function's return value is:
// =0: Intersection was not obtained.
// !=0: Intersection was obtained as follows:
//    =1: End point is a point on a perimeter of one of the surfaces.
//    =3: End point is one of boundary points in uvuv_list.
int MGPlane::isect_startHPL(
	const MGPosition& uvuv_startIn, //Starting point of the intersection line.
	MGPosition_list& uvuv_list,	//isect_startHPL will halt when ip reached one of 
		//the point in uvuv_list. isect_startHPL does not change uvuv_list(actually
		//uvuv_list is const). uvuv's space dimension is at least 4,
		//and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGSurface& sf2,	//2nd surface.
	MGSSisect& ssi,			//Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id
		//When the end point of ip was one of the points of uvuv_list, that is, 
		//when the function's return value was 3, uvuv_list's iterator
		//of the point will be returned,
		//When the end point was not a point of uvuv_list, end() of uvuv_list
		//will be returned.
) const{
		//cout<<"uvuv_startIn="<<uvuv_startIn<<" ,"<<uvuv_list<<endl;///////
	const MGPlane* pl2=dynamic_cast<const MGPlane*>(&sf2);
	if(!pl2){//Case that sf2 is not a plane.
	//Exchange this and surf parameter values in uvuv_list.
		MGPosition_list::iterator i=uvuv_list.begin(), ie=uvuv_list.end();
		for(;i!=ie; i++){
			double u=(*i)[0], v=(*i)[1];
			(*i).store_at(0,*i,2,2);
			(*i)(2)=u;(*i)(3)=v;
		}
	//Construct a working 1D surface sf1D2.
		MGSBRep* sf1D2=sf2.surf1D(*this);
	//Compute intersection lines using isect_startPlane.
		MGPosition uvuvS(uvuv_startIn);
		uvuvS(0)=uvuv_startIn[2]; uvuvS(1)=uvuv_startIn[3];
		uvuvS(2)=uvuv_startIn[0]; uvuvS(3)=uvuv_startIn[1];
		int obtained=sf2.isect_startPlane(uvuvS,uvuv_list,*this,ssi,uvuv_id);
		ssi.replace12();
		i=uvuv_list.begin(), ie=uvuv_list.end();
		for(;i!=ie; i++){
			double u=(*i)[0], v=(*i)[1];
			(*i).store_at(0,*i,2,2);
			(*i)(2)=u;(*i)(3)=v;
		}
		delete sf1D2;
		return obtained;
	}

	//Case that sf2 is also a plane.
	MGStraight sl;
	MGPSRELATION rl=relation(*pl2,sl);
	if(rl!=MGPSREL_ISECT)
		return 0;

	MGVector dir(sl.direction());
	if(uvuv_startIn.sdim()>4){
		MGVector dirin(3,uvuv_startIn,0,4);
		if(dirin%dir<0.)
			dir*=-1.;
	}
	MGPosition origin(eval(uvuv_startIn));
	sl.set_straight(MGSTRAIGHT_HALF_LIMIT, dir, origin);
	//Find the nearest point that lies on sl from uvuv_list,
	//which will be the end point.
	MGPosition_list::iterator i=uvuv_list.begin(), ie=uvuv_list.end();
	double tmin=-1.;
	for(; i!=ie; i++){
		double t;
		sl.on(eval(*i),t);
		if(tmin<0. || t<tmin){
			uvuv_id=i; tmin=t;
		}
	}
	if(tmin<0.){//If no end points found.
		MGPosition psuv1(2,uvuv_startIn), psuv2(2,uvuv_startIn,0,2);
		MGPosition pe=origin+dir;
		MGPosition uv;
		on(pe,uv); MGVector dir1(uv,psuv1);
		pl2->on(pe,uv); MGVector dir2(uv,psuv2);
		MGStraight sluv1(MGSTRAIGHT_HALF_LIMIT,dir1,psuv1);
		MGStraight sluv2(MGSTRAIGHT_HALF_LIMIT,dir2,psuv2);
		ssi=MGSSisect(sl,sluv1,sluv2);
		uvuv_id=ie;
		return 1;
	}

	//If an end point found.
	MGPosition pe(eval(*uvuv_id)), peuv1(2,*uvuv_id), peuv2(2,*uvuv_id,0,2);
	MGPosition psuv1(2,uvuv_startIn), psuv2(2,uvuv_startIn,0,2);
	ssi=MGSSisect(MGStraight(pe,origin),
						MGStraight(peuv1,psuv1),MGStraight(peuv2,psuv2));
	return 3;
}

//Intersection of Surface and a straight line.
MGCSisect_list MGPlane::isectSl(
	const MGStraight& sl,
	const MGBox& uvbox //indicates if this surface is restrictied to the parameter
					//range of uvbox. If uvbox.is_null(), no restriction.
)const{
	MGCSisect_list list(&sl,this);
	const MGBox& sbx=box();
	if(!sbx.crossing(sl))
		return list;

	MGCSisect is;
	if(sl.relation(*this,is)==MGPSREL_ISECT){
		const MGPosition& uv=is.param_surface();
		if(!uvbox.is_null()){
			if(!uvbox[0].includes(uv[0])) return list;
			if(!uvbox[1].includes(uv[1])) return list;
		}
		list.append(is);
	}
	return list;
}

// "isect_guess" computes one intersection point of surface and a curve,
// given initail guess parameter values of surface and curve.
//Function's return value is 1(true) if i.p. obtained, and 0(false) if not.
int MGPlane::isect_guess(
	const MGCurve& crv,		//Curve
	const MGPosition& uvi,	//Input initial guess parameter value
						// of the i.p. of the surface. 
	double ti,			//Input initial guess parameter value of the line.
	MGPosition& uv,		// Output parameter value obtained. 
	double& t)			// Output parameter value obtained. 
const{
	const MGStraight* sl=dynamic_cast<const MGStraight*>(&crv);
	if(sl) return isect_guess_straight(*sl,ti,uvi,t,uv);
	const MGCompositeCurve* ccrv=dynamic_cast<const MGCompositeCurve*>(&crv);
	if(ccrv) return isect_guess_composite(*ccrv,uvi,ti,uv,t);

	MGPosition A,P;
	MGVector AP,Lt,PQ;
	double dt,Lt2,AP2;
	MGStraight line; MGCSisect is;
	double error2=MGTolerance::wc_zero_sqr()*.5;
	double t0=crv.param_s(), t1=crv.param_e();

	int loop=0, tlow=0, thigh=0;
	t=ti; uv=uvi;
	while(loop++<16 && tlow<5 && thigh<5){
		A=eval(uv); P=crv.eval(t); AP=P-A;
	//A is guess point of the surface. P is guess point of the line.
		AP2=AP%AP;
		if(AP2<=error2)return 1;

		Lt=crv.eval(t,1);
		Lt2=Lt%Lt; if(MGMZero(Lt2)) break;
		line=MGStraight(MGSTRAIGHT_UNLIMIT,Lt,P);
		if(line.relation(*this,is)!=MGPSREL_ISECT) break;
		const MGPosition& Q=is.point();	//Q is the intersection point of osculating plane
	// of the surface at uv(A) and tangent straight line of the line at t(P).
		PQ=Q-P;
		uv=is.param_surface();
		dt=(PQ%Lt)/Lt2;
		
	// Update uv and t.
		t+=dt;// cout<<"guess dt="<<dt<<endl;//	
		if(t<t0)     {tlow+=1; thigh=0; t=t0;}
		else if(t>t1){tlow=0; thigh+=1; t=t1;}
		else         {tlow=0; thigh=0;}
	}	
	if(AP2<=error2*5.) return 1;
	return 0;
}

//与えられた曲線を自身に面直またはベクトル投影して曲線リストを求める。
//投影曲線は面上のパラメータ曲線と3次元曲線としてそれぞれ順番に、
//vec_crv_uv, vec_crvに格納される。
//uv曲線のトレランスはrc_zero()を、3次元曲線はline_zero()をそれぞれ使用している。
//get perpendicular or vector projection curve list.
//uv projection curves are put into vec_crv_uv(rc_zero() is used),
//3d projection curves are put into vec_crv(line_zero() is used) respectively.
//引数：
//		const MGCurve&			crv,		(I/ )	given curve.
//		MGPvector<MGCurve>&	vec_crv_uv,		( /O)	uv projection curve.
//		MGPvector<MGCurve>&	vec_crv,		( /O)	3d projection curve.
//		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
//戻り値：
//		投影曲線の数:		投影曲線が求まった
//		0:			投影曲線が求まらなかった
//		-1:			内部処理エラー
//		-2:			収束処理エラー（収束しなかった）
//追記：引数vecが与えられないとき、面直投影する。
//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int MGPlane::project(
	const MGCurve& crv,	//given curve.
	MGPvector<MGCurve>& vec_crv_uv,	//uv projection curve.
	MGPvector<MGCurve>& vec_crv,	//3d projection curve.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	assert(crv.sdim() > 1);
	const MGStraight* sl=dynamic_cast<const MGStraight*>(&crv);
	if(sl){
		return project(*sl,vec_crv_uv,vec_crv,vec);
	}
try{
	//初期化
	vec_crv_uv.clear(); vec_crv.clear();
	std::auto_ptr<MGCurve> crvRemoveKnot(crv.clone());
	crvRemoveKnot->remove_knot();

	if(vec == mgNULL_VEC ){	//面直投影
		//planeのノーマルでベクトル投影する
		return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, normal());
	}
	return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, vec);
}
catch(int err){
	return err;
}
}
int MGPlane::project(
	const MGStraight& sl,	//given curve.
	MGPvector<MGCurve>& vec_crv_uv,	//uv projection curve.
	MGPvector<MGCurve>& vec_crv,	//3d projection curve.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	vec_crv_uv.clear(); vec_crv.clear();
	MGSTRAIGHT_TYPE kind=sl.straight_type();
	if(kind==MGSTRAIGHT_EMPTY) return 0;

	const MGPosition P=sl.start_point();
	if(vec.is_null()){
		MGPosition uvP; on(P,uvP);
		if(kind==MGSTRAIGHT_SEGMENT){
			const MGPosition Q=sl.end_point();
			MGPosition uvQ; on(Q,uvQ);
			vec_crv_uv.push_back(new MGStraight(uvQ,uvP));
			vec_crv.push_back(new MGStraight(eval(uvQ),eval(uvP)));
			return 1;
		}
		const MGVector& dir=sl.direction();
		MGPosition uvdir; on(dir,uvdir);
		vec_crv_uv.push_back(new MGStraight(kind,uvdir,uvP));
		vec_crv.push_back(new MGStraight(kind,eval(uvdir),eval(uvP)));
	}else{
		MGStraight sl1(MGSTRAIGHT_UNLIMIT,vec,P);
		MGCSisect is1; relation(sl1,is1);
		const MGPosition& uvPprj=is1.param_surface();
		const MGPosition& Pprj=is1.point();
		if(kind==MGSTRAIGHT_SEGMENT){
			const MGPosition Q=sl.end_point();
			MGStraight sl2(MGSTRAIGHT_UNLIMIT,vec,Q);
			MGCSisect is2; relation(sl2,is2);
			const MGPosition& uvQprj=is2.param_surface();
			const MGPosition& Qprj=is2.point();
			vec_crv_uv.push_back(new MGStraight(uvQprj,uvPprj));
			vec_crv.push_back(new MGStraight(Qprj,Pprj));
			return 1;
		}
		const MGVector& dir=sl.direction();
		MGStraight sl3(MGSTRAIGHT_UNLIMIT,vec,dir);
		MGCSisect is3; relation(sl3,is3);
		const MGPosition& uvDirprj=is3.param_surface();
		const MGPosition& Dirprj=is3.point();
		vec_crv_uv.push_back(new MGStraight(kind,uvDirprj,uvPprj));
		vec_crv.push_back(new MGStraight(kind,Dirprj,Dirprj));
	}
	return 1;
}

//与えられた曲線を自身に面直またはベクトル投影して曲線リストを求める。
//投影曲線は3次元曲線としてvec_crvに格納される。
//3次元曲線のtoleranceはline_zero()を使用している。
//get perpendicular or vector projection curve list.
//3d projection curves are put into vec_crv(line_zero() is used).
//引数：
//		const MGCurve&			crv,		(I/ )	given curve.
//		MGPvector<MGCurve>&	vec_crv,		( /O)	3d projection curve.
//		const MGVector&			vec=mgNULL_VEC	(I/ )	projection vector.
//戻り値：
//		投影曲線の数:		投影曲線が求まった
//		0:			投影曲線が求まらなかった
//		-1:			内部処理エラー
//		-2:			収束処理エラー（収束しなかった）
//追記：引数vecが与えられないとき、面直投影する。
//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is general world coordinate curves('vec_crv')
int MGPlane::project(
	const MGCurve& crv,	//given curve.
	MGPvector<MGCurve>& vec_crv,//3d projection curve.
	const MGVector& vec	//projection vector.
						//if vec = NULL then calculate perpendicular project.
)const{
	MGPvector<MGCurve> vec_crv_uv;
	return project(crv,vec_crv_uv,vec_crv,vec);
}
