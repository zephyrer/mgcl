/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/KnotVector.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/SurfCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/BSumCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/RLBRep.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGTrimmedCurve Class.
//MGTrimmedCurve is a part of original curve that has limitted parameter
//range.
//MGTrimmedCurve is a temporal curve, and does not have update functions.
//

///////Constructor///////

MGTrimmedCurve::MGTrimmedCurve()
:MGCurve(),m_sameRange(false),m_curve(0), m_knotV(0){;}

//Copy constructor.
MGTrimmedCurve::MGTrimmedCurve(const MGTrimmedCurve& tc)
:MGCurve(tc), m_knotV(0), m_curve(tc.m_curve),
m_sameRange(tc.m_sameRange),m_range(tc.m_range){
	if(tc.m_knotV)
		m_knotV=new MGKnotVector(*(tc.m_knotV));
}

MGTrimmedCurve::MGTrimmedCurve(const MGCurve& crv, double t1, double t2)
:m_curve(&crv),m_sameRange(false),m_range(t1,t2), m_knotV(0){
	update_mark();
	if(t1>t2)
		m_range=MGInterval(t2,t1);
	const MGInterval& crange=crv.param_range();
	m_range&=crange;
	if(m_range==crange){
		m_range=crange;
		m_sameRange=true;
	}
	const MGTrimmedCurve* tc=dynamic_cast<const MGTrimmedCurve*>(&crv);
	if(tc) m_curve=tc->m_curve;
}

MGTrimmedCurve::MGTrimmedCurve(const MGCurve& crv, const MGInterval range)
:m_curve(&crv),m_sameRange(false),m_range(range),m_knotV(0){
	update_mark();
	const MGInterval& crange=crv.param_range();
	m_range&=crange;
	if(m_range==crange){
		m_range=crange;
		m_sameRange=true;
	}
	const MGTrimmedCurve* tc=dynamic_cast<const MGTrimmedCurve*>(&crv);
	if(tc) m_curve=tc->m_curve;
}

//////////Destructor//////////
MGTrimmedCurve::~MGTrimmedCurve(){
	if(m_knotV) delete m_knotV;
}

//Assignment.
//When the leaf object of this and obj2 are not equal, this assignment
//does nothing.
MGTrimmedCurve& MGTrimmedCurve::operator=(const MGTrimmedCurve& tc){
	if(this==&tc)
		return *this;

	MGCurve::operator=(tc);
	if(m_knotV){
		delete m_knotV;
		m_knotV=0;
	}
	if(tc.m_knotV)
		m_knotV=new MGKnotVector(*(tc.m_knotV));
	m_curve=tc.m_curve;
	m_sameRange=tc.m_sameRange;
	m_range=tc.m_range;
	return *this;
}
MGTrimmedCurve& MGTrimmedCurve::operator=(const MGGel& gel2){
	const MGTrimmedCurve* gel2_is_this=dynamic_cast<const MGTrimmedCurve*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

///////Member Function///////

//Returns B-Rep Dimension.
size_t MGTrimmedCurve::bdim() const{
	const MGKnotVector& t=knot_vector();
	return t.bdim();
}

//Return minimum box that includes the curve of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox MGTrimmedCurve::box_limitted(
	const MGInterval& rng // Parameter Range of the curve.
)const{
	MGInterval prng(param_range()&rng);
	return m_curve->box_limitted(prng);
}

//Compute the box of the whole of the curve.
//Returned is a newed object pointer.
MGBox* MGTrimmedCurve::compute_box()const{
	MGBox cbox=m_curve->box_limitted(m_range);
	return new MGBox(cbox);
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGCurve* MGTrimmedCurve::clone() const{
	if(m_sameRange)
		return m_curve->clone();
	return m_curve->part(m_range.low_point(),m_range.high_point());
}

//Exclusive function for common.
//When this curve is a TrimmedCurve of MGComposite, narrow the parameter range
//by this m_range.
void MGTrimmedCurve::narrow_into_range(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
)const{
	std::vector<double> paramVec;
	size_t n=vecComSpan.size();
	for(size_t i=0; i<n;){
		double& s0=vecComSpan[i++];
		double& s1=vecComSpan[i++];
		double& t0=vecComSpan[i++];
		double& t1=vecComSpan[i++];
		if(s1<=m_range)
			continue;
		if(m_range<=s0)
			continue;
		if(m_range>>s0 && m_range>>s1){//when span [s0, s1] is included in m_range
			paramVec.push_back(s0);
			paramVec.push_back(s1);
			paramVec.push_back(t0);
			paramVec.push_back(t1);
			continue;
		}

		double s0d=m_range.low_point();
		double s1d=m_range.high_point();
		double t0d=t0, t1d=t1;
		if(s0<=s0d && s0d<=s1){
			double tg=t0+(t1-t0)*((s0d-s0)/(s1-s0));
			curve2.perp_guess(t0,t1,start_point(),tg,t0d);
		}
		if(s0<=s1d && s1d<=s1){
			double tg=t0+(t1-t0)*((s1d-s0)/(s1-s0));
			curve2.perp_guess(t0,t1,end_point(),tg,t1d);
		}
		paramVec.push_back(s0d);
		paramVec.push_back(s1d);
		paramVec.push_back(t0d);
		paramVec.push_back(t1d);
	}
	vecComSpan=paramVec;
}

//関数名：common
//目的：与えられた曲線と自身の共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			curve2,		(I/ )	与えられる曲線
//		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
//		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
//					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
//				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
//		MGCCisect_list&			isect		( /O)	交点
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
int MGTrimmedCurve::common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan,
	MGCCisect_list& isect
)const{
	int ret;
	const MGCompositeCurve* cmpo=dynamic_cast<const MGCompositeCurve*>(m_curve);
	if(cmpo){
		ret=cmpo->common(curve2,vecComSpan,isect);
	}else{
		ret=m_curve->common(curve2,vecComSpan,isect);
	}
	if(!m_sameRange)
		narrow_into_range(curve2,vecComSpan);

	MGCCisect_list::CCiterator j=isect.begin(), jend=isect.end(),j2;
	while(j!=jend){//Remove parameter values that are outside the range.
		j2=j; j2++;
		if(m_range<<(*j).param1()) isect.removeAt(j);
		j=j2;
	}
	return ret;
}

//関数名：common
//目的：与えられた曲線と自身の共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			curve2,		(I/ )	与えられる曲線
//		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
//		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
//					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
//				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
//戻り値：
//		共通部分の数:	共通部分が求まった
//		0:				共通部分がなかった
//		-1:				共通エッジの収束計算エラー
//		-2:				共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
//  パラメータ範囲*rc_zero()を使用した
int MGTrimmedCurve::common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
)const{
	int ret;
	const MGCompositeCurve* cmpo=dynamic_cast<const MGCompositeCurve*>(m_curve);
	if(cmpo){
		ret=cmpo->common(curve2,vecComSpan);
	}else{
		ret=m_curve->common(curve2,vecComSpan);
	}
	if(!m_sameRange)
		narrow_into_range(curve2,vecComSpan);
	return ret;
}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGTrimmedCurve::copy_as_nurbs() const{
	MGCurve* tcrv=m_curve->part(m_range.low_point(), m_range.high_point());
	MGCurve* nurbs=tcrv->copy_as_nurbs();
	delete tcrv;
	return nurbs;
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGCurve* MGTrimmedCurve::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2 		// Source order of this line.
)const{
	MGCurve* crv=(MGCurve*)
		(m_curve->copy_change_dimension(sdim,start1,start2));
	MGCurve* trimmed=crv->part(m_range.low_point(), m_range.high_point());
	delete crv;
	return trimmed;
}

//Construct new curve object by copying to newed area,
//and limitting the parameter range to prange.
//Returned is newed object and must be deleted.
MGCurve* MGTrimmedCurve::copy_limitted(const MGInterval& prange) const{
	MGInterval prng2(prange&m_range);
	return m_curve->part(prng2.low_point(), prng2.high_point());
}
MGCurve* MGTrimmedCurve::copy_limitted() const{
	return m_curve->part(m_range.low_point(), m_range.high_point());
}

//Compute curvilinear integral of the 1st two coordinates.
//This integral can be used to compute area sorounded by the curve.
//(線積分）を求める。
//curvilinear_integral from t1 to t2 can be obtained by
//Integral of (x*dy-y*dx) about t, where curve is expressed by
//f(t)=(x(t),y(t)), dx=dx/dt, and dy=dy/dt.
double MGTrimmedCurve::curvilinear_integral(double t1, double t2) const{
	MGInterval prng2(m_range&MGInterval(t1,t2));
	return m_curve->curvilinear_integral
					(prng2.low_point(), prng2.high_point());
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGTrimmedCurve::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	double lenPara;
	double t, ts=m_curve->param_s(), te=m_curve->param_e();
	if(start){
		t=param_s();
		lenPara=m_curve->length(ts,t);
	}else{
		t=param_e();
		lenPara=m_curve->length(t,te);
	}
	if(length>=lenPara)
		return;

	if(start)
		length*=-1.;
	double t2=m_curve->length_param(t,length);
	if(start)
		change_range(t2,param_e());
	else
		change_range(param_s(),t2);
}

//Provide divide number of curve span for function intersect.
size_t MGTrimmedCurve::intersect_dnum() const{
	size_t n;
	const MGKnotVector& t=knot_vector();
	if(t==mgNULL_KNOT_VECTOR) n=m_curve->intersect_dnum();
	else{
		size_t k=order();
		int nspan=bdim()+2-k;
		int km2=k-2; if(km2<=0) km2=1;
		double t1=m_range.low_point(), t2=m_range.high_point();
		size_t i1=t.locate(t1), i2=t.locate(t2,1);
		if(i2<=i1)
			i2=i1+1;
		if(dynamic_cast<const MGRLBRep*>(m_curve)) n=(i2-i1+1)*(km2+1);
		else n=(i2-i1+1)*km2;
	}
	return n;
}

//Intersection of Curve.
MGCCisect_list MGTrimmedCurve::isect_in_range(const MGCurve& crv2) const{
	MGCCisect_list list=m_curve->isect(crv2);
	if(m_sameRange) return list;
	MGCCisect_list::CCiterator i=list.begin(), iend=list.end(),i2;
	while(i!=iend){//Remove parameter values that are outside the range.
		i2=i; i2++;
		if(m_range<<(*i).param1()) list.removeAt(i);
		i=i2;
	}
	return list;
}

//Intersection of Curve and Surface
MGCCisect_list MGTrimmedCurve::isect(const MGCurve& curve2) const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGStraight& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGRLBRep& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGEllipse& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGLBRep& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGSurfCurve& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGBSumCurve& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGTrimmedCurve& curve2)const{
	return isect_in_range(curve2);
}

MGCCisect_list MGTrimmedCurve::isect(const MGCompositeCurve& curve2)const{
	return isect_in_range(curve2);
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list MGTrimmedCurve::isect_withC1LB(const MGLBRep& curve2)const{
	return isect_in_range(curve2);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGSurface& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGPlane& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGSphere& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGCylinder& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGSBRep& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGRSBRep& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect(const MGBSumSurf& surf) const{
	return isect_in_range(surf);
}

//Intersection with a surface.
MGCSisect_list MGTrimmedCurve::isect_in_range(const MGSurface& surf)const{
	MGCSisect_list list=m_curve->isect(surf);
	if(m_sameRange)
		return list;
	MGCSisect_list::CSiterator i=list.begin(), iend=list.end(),i2;
	while(i!=iend){//Remove parameter values that are outside the range.
		i2=i; i2++;
		if(m_range<<(*i).param_curve())
			list.removeAt(i);
		i=i2;
	}
	return list;
}

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
MGCParam_list MGTrimmedCurve::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
)const{
	MGCParam_list list=m_curve->isect_1D(f,coordinate);
	if(m_sameRange) return list;
	MGCParam_list::Citerator i=list.begin(), iend=list.end(),i2;
	while(i!=iend){//Remove parameter values that are outside the range.
		i2=i; i2++;
		if(m_range<<(*i)) list.removeAt(i);
		i=i2;
	}
	return list;
}

//Cmpute curve length of the interval.
//If t1 is greater than t2, return negative value.
// 与えられたパラメータ値間の曲線の長さを返す。
// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
double MGTrimmedCurve::length(double t1, double t2) const{
	t1=range(t1); t2=range(t2);
	return m_curve->length(t1,t2);
}

//Inverse function of length. Compute the point that is away from
//the point t by length len.
// lengthの逆関数。指定パラメータtで示される点から指定距離len
// 曲線上に沿って離れた点を示すパラメータ値を返す。
double MGTrimmedCurve::length_param( double t, double len) const{
	t=range(t);
	return m_curve->length_param(t,len);
}

//Update this by limiting the parameter range of the curve.
// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
MGTrimmedCurve& MGTrimmedCurve::limit(const MGInterval& nrng){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	m_range&=nrng;
	if(m_range==m_curve->param_range()){
		m_range=m_curve->param_range();
		m_sameRange=true;
	}else m_sameRange=false;
	update_mark();
	return *this;
}

// Return ending parameter value.
double MGTrimmedCurve::param_e() const{
	return m_range.high_point();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGTrimmedCurve::param_normalize(double t) const{
	double s1=m_range.low_point(); if(t<=s1) return s1;
	double s2=m_range.high_point(); if(t>=s2) return s2;
	double len=s2-s1;
	if((t-s1)<=(s2-t)){if(MGREqual_base(t,s1,len)) return s1;}
	else {if(MGREqual_base(t,s2,len)) return s2;}
	return m_curve->param_normalize(t);
}

// Return starting parameter value.
double MGTrimmedCurve::param_s() const{
	return m_range.low_point();
}
	
//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGCurve* MGTrimmedCurve::part(double t1, double t2, int multiple) const{
	assert(t1<t2);
	MGInterval rng(param_range()&MGInterval(t1,t2));
	return m_curve->part(rng.low_point(),rng.high_point(),multiple);
}

//Compute all foot points of the perpendicular line from point to
//the curve.
// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
// すべて求める。
MGCParam_list MGTrimmedCurve::perps(
	const MGPosition& P		//Point(指定点)
)const{
	if(m_sameRange) return m_curve->perps(P);
	MGCParam_list list=m_curve->perps(P);
	MGCParam_list::Citerator i=list.begin(), iend=list.end(),i2;
	while(i!=iend){//Remove parameter values that are outside the range.
		i2=i; i2++;
		if(m_range<<(*i)) list.removeAt(i);
		i=i2;
	}
	return list;
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGTrimmedCurve::perps(const MGSurfCurve& crv2)const{
	return perps_in_range(crv2);
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGTrimmedCurve::perps_in_range(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list=m_curve->perps(crv2);
	if(m_sameRange)
		return list;

	MGPosition_list::iterator i=list.begin(), iend=list.end(),i2;
	while(i!=iend){//Remove parameter values that are outside the range.
		i2=i; i2++;
		if(m_range<<((*i).ref(0))) list.removeAt(i);
		i=i2;
	}
	return list;
}

//Round t into curve's parameter range.
// 入力パラメータをパラメータ範囲でまるめて返却する。
double MGTrimmedCurve::range(double t) const{
	double s=m_range.low_point();
	if(t<=s) return s;
	s=m_range.high_point();
	if(t>=s) return s;
	return t;
}

//Unlimit parameter range of the curve(limitをはずす)
MGCurve& MGTrimmedCurve::unlimit(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	m_range=m_curve->param_range();
	m_sameRange=true;
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGTrimmedCurve::unlimit_end(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	m_range(1)=m_curve->param_range()[1];
	if(m_range==m_curve->param_range()){
		m_range=m_curve->param_range();
		m_sameRange=true;
	}else m_sameRange=false;
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGTrimmedCurve::unlimit_start(){
	if(m_knotV){
		delete m_knotV; m_knotV=0;
	}

	m_range(0)=m_curve->param_range()[0];
	if(m_range==m_curve->param_range()){
		m_range=m_curve->param_range();
		m_sameRange=true;
	}else
		m_sameRange=false;
	update_mark();
	return *this;
}

//Operator overload(演算子多重定義)

//Logical operator overload(論理演算子多重定義)

// 与曲線と自身が等しいかの比較判定を行う。
bool MGTrimmedCurve::operator==(const MGTrimmedCurve& crv2) const{
	if((crv2.m_range)!=m_range)
		return false;
	if(m_curve==0){
		if(crv2.m_curve==0)
			return true;
		return false;
	}
	if(crv2.m_curve==0){
		if(m_curve==0)
			return true;
		return false;
	}

	return *(crv2.m_curve)==(*m_curve);
}
bool MGTrimmedCurve::operator==(const MGCompositeCurve& crv2) const{
	return is_same_curve(crv2);
}
bool MGTrimmedCurve::operator<(const MGTrimmedCurve& gel2)const{
	if(*(gel2.m_curve)==*m_curve)
		return m_range.length()<gel2.m_range.length();
	return *m_curve<*(gel2.m_curve);
}
bool MGTrimmedCurve::operator==(const MGGel& gel2)const{
	const MGTrimmedCurve* gel2_is_this=dynamic_cast<const MGTrimmedCurve*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	else{
		const MGCompositeCurve* gel2_is_compo=dynamic_cast<const MGCompositeCurve*>(&gel2);
		if(gel2_is_compo)
			return operator==(*gel2_is_compo);

	}
	return false;
}
bool MGTrimmedCurve::operator<(const MGGel& gel2)const{
	const MGTrimmedCurve* gel2_is_this=dynamic_cast<const MGTrimmedCurve*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

bool MGTrimmedCurve::is_same_curve(const MGCurve& crv2)const{
	if(!m_sameRange)
		return false;
	if(m_curve==0)
		return false;
	return (*m_curve)==(crv2);
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGTrimmedCurve::is_planar(MGPlane& plane)const{
	if(m_sameRange)
		return m_curve->is_planar(plane);
	MGCurve* tcrv=this->part(m_range.low_point(), m_range.high_point());
	bool planar=tcrv->is_planar(plane);
	delete tcrv;
	return planar;
}

//Access to i-th element of knot
double MGTrimmedCurve::knot(size_t i) const{
	const MGKnotVector& t=knot_vector();
	return t(i);
}

//Returns the knot vector.
const MGKnotVector& MGTrimmedCurve::knot_vector() const{
	if(m_sameRange) return m_curve->knot_vector();
	if(!m_knotV){
		m_knotV=new MGKnotVector(m_curve->knot_vector(),
					m_range.low_point(), m_range.high_point());
	}
	return *m_knotV;
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGTrimmedCurve& MGTrimmedCurve::coordinate_exchange(size_t i, size_t j){
	assert(false); return *this;
}

//Negate the curve direction(曲線の方向を反転する)
void MGTrimmedCurve::negate(){
	assert(false);
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGTrimmedCurve::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
) const{
	std::auto_ptr<MGCurve> crv=m_curve->oneD(g);
	crv->limit(m_range);
	return crv;
}

//ノット削除関数(B表現曲線のみ)
//トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
//戻り値は、削除したノットの数
//removal knot. line_zero tolerance is used.
void MGTrimmedCurve::remove_knot(){assert(false);}

//Update the curve by translation.
// 与ベクトルだけ曲線を平行移動して自身とする。
MGTrimmedCurve& MGTrimmedCurve::operator+= (const MGVector& v)
{assert(false); return *this;}

//Update the curve by translation.
// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGTrimmedCurve& MGTrimmedCurve::operator-= (const MGVector& v)
{assert(false); return *this;}

//Update the curve by multiplying scale.
// 与えられたスケールを曲線にかける。
MGTrimmedCurve& MGTrimmedCurve::operator*= (double scale)
{assert(false); return *this;}

//Update the curve by transformation of matrix.
// 与えられた変換で直線の変換を行い自身の直線とする。
MGTrimmedCurve& MGTrimmedCurve::operator*= (const MGMatrix& mat)
{assert(false); return *this;}

//Update the curve by transformation of transf.
// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGTrimmedCurve& MGTrimmedCurve::operator*= (const MGTransf& tr)
{assert(false); return *this;}

MGTrimmedCurve MGTrimmedCurve::operator+ (const MGVector& v) const{
	assert(false); return *this;
}
MGTrimmedCurve operator+ (const MGVector& v, const MGTrimmedCurve& cv2){
	assert(false); return cv2;
}
MGTrimmedCurve MGTrimmedCurve::operator- (const MGVector& v) const{
	assert(false); return *this;
}
MGTrimmedCurve MGTrimmedCurve::operator* (double scale) const{
	assert(false); return *this;
}
MGTrimmedCurve operator* (double scale, const MGTrimmedCurve& cv2){
	assert(false); return cv2;
}
MGTrimmedCurve MGTrimmedCurve::operator* (const MGMatrix& mat) const{
	assert(false); return *this;
}
MGTrimmedCurve MGTrimmedCurve::operator* (const MGTransf& tr) const{
	assert(false); return *this;
}
