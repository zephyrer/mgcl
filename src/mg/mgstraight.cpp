/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Transf.h"
#include "mg/Unit_vector.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/CSisect.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Ellipse.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGStraighy Class Implementation
//

//
// Constructor
	
//Void constructor
MGStraight::MGStraight()
: MGCurve()	,m_sparam(0.0),m_endparam(-1.),m_knotV(0){;}

//Straight specifying all the member data. All of the data are employed as 
//member data of this straight.
MGStraight::MGStraight (
	const MGEReal& endparam,	//end parameter value
	const MGEReal& sparam,		//start parameter value
	const MGVector& direction,	//Direction, which will be the direction of this.
	const MGPosition& Origin	//Origin
):m_direction(direction), m_endparam(endparam),
m_knotV(0), m_root_point(Origin), m_sparam(sparam){;}

//Straight specifying all the member data. All of the data are employed as 
//member data of this straight.
MGStraight::MGStraight (
	double endparam,			//end parameter value
	double sparam,				//start parameter value
	const MGVector& direction,	//Direction, which will be the direction of this.
	const MGPosition& Origin	//Origin
):m_direction(direction), m_endparam(endparam),
m_knotV(0), m_root_point(Origin), m_sparam(sparam){;}

//Copy constructor.
MGStraight::MGStraight(const MGStraight& sl)
:MGCurve(sl), m_direction(sl.m_direction), m_endparam(sl.m_endparam),
m_knotV(0), m_root_point(sl.m_root_point), m_sparam(sl.m_sparam){;}

// 始点、方向ベクトル、直線のタイプを指定して直線を生成する。
MGStraight::MGStraight (
	MGSTRAIGHT_TYPE t,
	const MGVector& v,
	const MGPosition& p
):MGCurve(),m_root_point(p),m_direction(v.normalize())
,m_sparam(0.0),m_endparam(MGINFINITE_PLUS),m_knotV(0){
	if (t==MGSTRAIGHT_SEGMENT) m_endparam=v.len();
	else if(t==MGSTRAIGHT_EMPTY) m_endparam=-1.;
	else if(t==MGSTRAIGHT_UNLIMIT) m_sparam=MGEReal(MGINFINITE_MINUS);

	size_t dimv=v.sdim(), dimp=p.sdim();
	if(dimv<dimp) m_direction=MGVector(dimp,v);
	else if(dimv>dimp) m_root_point=MGVector(dimv,p);
}

// ２点から直線を生成する。
//Parameter value of the start point is set to be 0.
MGStraight::MGStraight(
	const MGPosition& end,		//End point.
	const MGPosition& start		//Start point.
):MGCurve(),m_root_point(start),m_sparam(0.0),m_knotV(0){
	// ２点から単位ベクトルを作成する。
	MGVector v(end,start);
	m_direction = v.normalize();

	// 単位方向ベクトルをもつので
	// 作成したベクトルの長さが終点のパラメータ値
	m_endparam=v.len();
	if(m_root_point.sdim()<m_direction.sdim())
		m_root_point=MGVector(m_direction.sdim(),start);
}

//始終点の座標、パラメータ値から直線を生成する。
//MGSTRAIGHT_SEGMENT straight from two points.
//Start point is start and end point is end.
//In this version, can specify start and end parameter values.
MGStraight::MGStraight (
	const MGPosition& Pe,const MGPosition& Ps, 
	const double Te,	const double Ts
):MGCurve(),m_sparam(Ts),m_endparam(Te),m_knotV(0){
	assert(Ts<Te);
	// ２点から方向ベクトルを作成する。
	MGVector v(Pe,Ps);
	m_direction = v/(Te-Ts);
	m_root_point = Ps - Ts * m_direction;
}

// 始点、単位方向ベクトル、終点のパラメータ値を指定して直線を生成する。
MGStraight::MGStraight (
	const MGUnit_vector & v,
	double d,
	const MGPosition& p
):MGCurve(),m_direction(v),m_sparam(0.0),m_endparam(d),
	m_root_point(p),m_knotV(0){
	size_t dimv=v.sdim(), dimp=p.sdim();
	if(dimv<dimp) m_direction=MGVector(dimp,v);
	else if(dimv>dimp) m_root_point=MGVector(dimv,p);
	if (m_endparam < 0.0) {
		m_direction = -m_direction;
		m_endparam = -m_endparam;
	}
}
	
//Construct the infinite straight line that is a perpendicular bisect
//of the two point P1 and P2 and that is normal to the vector N.
//The line's direction is N*(P2-P1).
//N is the normal of the plane P1, P2, and the constructed line lie on.
MGStraight::MGStraight (
	const MGPosition& P1,	//point 1.
	const MGPosition& P2,	//point 2.
	const MGVector& N
):MGCurve(),m_root_point((P1+P2)*.5),m_direction(N*(P2-P1))
,m_sparam(MGINFINITE_MINUS),m_endparam(MGINFINITE_PLUS),m_knotV(0){
}

//Construct Straight Line copying original line. Able to change
//space dimension and ordering of axis.
MGStraight::MGStraight (
	size_t dim,				  //New space dimension.
	const MGStraight& line2,  //Original line.
	size_t start1, 			  //store order of new line.
	size_t start2)	 		  //source order of original line.
:MGCurve(line2), m_root_point(dim,line2.m_root_point,start1, start2)
	, m_direction(dim, line2.m_direction, start1,start2)
	,m_sparam(line2.m_sparam), m_endparam(line2.m_endparam),m_knotV(0){
	update_mark();
	save_length_zero();	
}

///////// Destructor ///////////
MGStraight::~MGStraight(){
	if(m_knotV) delete m_knotV;
}

//
// メンバ関数
//
// 指定線分を囲むボックスを返却する。
MGBox MGStraight::box_limitted(const MGInterval& l) const{
	// 直線範囲のパラメータを作成して入力Intervalとの積をとる。
	MGStraight sl(*this); sl.limit(l);
	return sl.box();
}

//Changing this object's space dimension.
MGStraight& MGStraight::change_dimension(
	size_t dim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	m_root_point=MGPosition(dim,m_root_point,start1, start2);
	m_direction=MGVector(dim, m_direction, start1,start2);
	save_length_zero();	
	update_mark();
	return *this;
}

//Change parameter range, be able to change the direction by providing
//s1 greater than s2.
void MGStraight::change_range(
	double s1,		//Parameter value for the start of the original. 
	double s2)		//Parameter value for the end of the original. 
{
	if(m_knotV){delete m_knotV; m_knotV=0;}
//	cout<<eval(param_s())<<","<<eval(param_e())<<endl;/////
	if(s1>s2){
		negate();
		double save=s1; s1=s2; s2=save;
	}
	if(m_sparam.minus_infinite()&&m_endparam.plus_infinite())
		return;

	double t1=m_sparam.value(), t2=m_endparam.value();
	if(m_sparam.minus_infinite()){
		m_root_point+=(t2-s2)*m_direction;
		m_endparam=s2;
		return;
	}
	if(m_endparam.plus_infinite()){
		m_root_point+=(t1-s1)*m_direction;
		m_sparam=s1;
		return;
	}
	double s2ms1=s2-s1, t2mt1=t2-t1;
	if(MGMZero(s2ms1)){
		if(MGMZero(t2mt1)){
			m_root_point+=(t2-s2)*m_direction;
			m_sparam=s1; m_endparam=s2;
		}
		return;
	}
	double ratio=t2mt1/s2ms1;
	m_root_point += m_direction*(t1-s1*ratio);
	m_direction *= ratio;
	m_sparam=s1; m_endparam=s2;
//	cout<<eval(param_s())<<","<<eval(param_e())<<endl;/////
}

//Compute box of the whole line.
MGBox* MGStraight::compute_box() const{
	MGInterval p_limit=param_range();
	if(p_limit.empty()) return new MGBox();

	// Intervalで示されるboxを作成する。
	if(p_limit.finite())
		return new MGBox(eval_position(p_limit.low_point()),
					eval_position(p_limit.high_point()) );

	size_t dim=sdim();
	MGBox* box=new MGBox(dim);
	if (p_limit.finite_below()&&p_limit.infinite_above()) {
		for(size_t i=0; i<dim; i++){
			if(MGRZero2(m_direction.ref(i),m_direction.len()))
				(*box)(i)=MGInterval(start_point().ref(i),start_point().ref(i));
			else if(m_direction.ref(i)>0.)
				(*box)(i)=
				MGInterval(MGINTERVAL_FINITE_BELOW,start_point().ref(i));
			else
				(*box)(i)=
				MGInterval(MGINTERVAL_FINITE_ABOVE,start_point().ref(i));
		}
	} else if (p_limit.finite_above()&&p_limit.infinite_below()) {
		for(size_t i=0; i<dim; i++){
			if(MGRZero2(m_direction.ref(i),m_direction.len()))
				(*box)(i)=MGInterval(root_point().ref(i),root_point().ref(i));
			else if(m_direction.ref(i)>0.)
				(*box)(i)=
				MGInterval(MGINTERVAL_FINITE_ABOVE,root_point().ref(i));
			else
				(*box)(i)=
				MGInterval(MGINTERVAL_FINITE_BELOW,root_point().ref(i));
		}
	} else {
		for(size_t i=0; i<dim; i++){
			if(MGRZero2(m_direction.ref(i),m_direction.len())){
				MGInterval Ii=MGInterval(root_point().ref(i),root_point().ref(i));
				(*box)(i)=Ii;
			} else{
				MGInterval Ii=MGInterval(MGINTERVAL_INFINITE);
				(*box)(i)=Ii;
			}
		}
	} 
	return box;
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGStraight& MGStraight::coordinate_exchange(size_t i, size_t j){
	assert(i<sdim() && j<sdim());
	m_direction.swap(i,j); m_root_point.swap(i,j);
	update_mark();
	return *this;
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGStraight* MGStraight::clone() const{return new MGStraight(*this);}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGStraight::copy_as_nurbs() const{
	return new MGLBRep(*this);
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGStraight* MGStraight::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2 		// Source order of this line.
)const{
	return new MGStraight(sdim,*this,start1,start2);
}
//Compute curvilinear integral of the 1st two coordinates.
//This integral can be used to compute area sorounded by the curve.
//線積分を求める。
double MGStraight::curvilinear_integral(double t1, double t2) const{
	t1=range(t1); t2=range(t2);
	return (t2-t1)*(root_point().ref(0)*m_direction.ref(1)
				  - root_point().ref(1)*m_direction.ref(0));
}

// 自身と与えられた点との距離を返す。
double MGStraight::distance(
	const MGPosition& p ) const
{
	if(straight_type() == MGSTRAIGHT_EMPTY) return 0.0;

	MGVector v1(p,root_point());  // 基点から与点へのベクトル。
	double v1len=v1.len();
	double dlen=m_direction.len();
	double hs = (v1%m_direction)/dlen;
						// 与点から直線への垂線の足と基点との距離。
	double t=hs/dlen;	//parameter of the straight.

	double d;
	if(t<m_sparam) d=(p-start_point()).len();
	else if(t>m_endparam) d=(p-end_point()).len();
	else{
		d= v1len*v1len - hs*hs; if(d<0.) d=0.;
		d= sqrt(d);	// 点と無限直線の最短距離。
	}
	// 最短距離を返す。
	return d;
}

// 自身と与えられた直線との距離を返す。
double MGStraight::distance(const MGStraight& s) const{
	if(straight_type()==MGSTRAIGHT_EMPTY
		|| s.straight_type()==MGSTRAIGHT_EMPTY) return 0.;

	// 与えられた直線との関係を調べる。
	MGCCisect inter;
	MGPSRELATION rl=relation(s,inter);
	if(rl==MGPSREL_ISECT) return 0.;	// 交差のとき距離は０

	const MGStraight* sl1; const MGStraight* sl2;
	if(straight_type()==MGSTRAIGHT_SEGMENT)          {sl1=this; sl2=&s;}
	else if(s.straight_type()==MGSTRAIGHT_SEGMENT)   {sl2=this; sl1=&s;}
	else if(straight_type()==MGSTRAIGHT_HALF_LIMIT)  {sl1=this; sl2=&s;}
	else if(s.straight_type()==MGSTRAIGHT_HALF_LIMIT){sl2=this; sl1=&s;}
	else{	//When both are infinite line.
		if (rl==MGPSREL_TORSION){
			MGPlane pl(m_direction*s.m_direction, s.root_point());
			return pl.distance(root_point());
		}else{//When parallel or coincidence.
			return distance(s.root_point());
		}
	}

	MGPosition_list list=perps(s);
	if(list.size()){
		MGVector v1=eval((list.first()).ref(0));
		MGVector v2=s.eval((list.first()).ref(1));
		return (v1-v2).len();
	}

	double d1s=sl2->distance(sl1->start_point()),
		d1e=sl2->distance(sl1->end_point()),
		d2s=sl1->distance(sl2->start_point()),
		d2e=sl1->distance(sl2->end_point());

	double d;
	if(sl1->straight_type()==MGSTRAIGHT_SEGMENT){
		d=d1s; if(d>d1e) d=d1e;
		if(sl2->straight_type()==MGSTRAIGHT_SEGMENT){
			if(d>d2s) d=d2s; if(d>d2e) d=d2e;
		}else if(sl2->straight_type()==MGSTRAIGHT_HALF_LIMIT){
			if((sl2->m_sparam).finite()){ if(d>d2s) d=d2s;}
			else{if(d>d2e) d=d2e;}
		}
	}else if(sl2->straight_type()==MGSTRAIGHT_HALF_LIMIT){
		if((sl1->m_sparam).finite()){
			if((sl2->m_sparam).finite()){
				if(d1s<=d2s) d=d1s; else d=d2s;
			}else{
				if(d1s<=d2e) d=d1s; else d=d2e;
			}
		}else{
			if((sl2->m_sparam).finite()){
				if(d1e<=d2s) d=d1e; else d=d2s;
			}else{
				if(d1e<=d2e) d=d1e; else d=d2e;
			}
		}
	}else{//When sl1==half limit and sl2==ssegment.
		d=d2s; if(d>d2e) d=d2e;
	}
	// 最短距離を返す
	return d;
}

// 終点の座標値を返却する。
MGPosition MGStraight::end_point() const{
	return eval_position(m_endparam.value());
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGStraight::eval(
	double t,		// Parameter value.
	size_t nderiv,	// Order of Derivative.
	int left		//Left continuous(left=true)
					//or right continuous(left=false).
	) const{
	if(nderiv==0)      return eval_position(t);
	else if(nderiv==1) return m_direction;
	else               return MGVector(sdim(),MGVector (0.0,0.0));
}

// パラメータ値を与えて位置、一次微分値、二次微分値を求める。
void MGStraight::eval_all (double d,
	MGPosition & p,				   // 位置      
	MGVector & v1,				   // 一次微分値
	MGVector & v2				   // 二次微分値
	) const {
	p = eval_position (d);			// 位置
	v1 = m_direction;				// 一次微分値
	v2 = MGVector(sdim(),0.0);		// 二次微分値
}

// 直線上の与えられたパラメータ値における一次微分値を求める。
MGVector MGStraight::eval_deriv(double d) const{
	return m_direction;
}

// 与えられたパラメータ値に相当する直線上の点を返却する。
MGPosition MGStraight::eval_position (double d)const{
	return m_root_point+range(d)*m_direction;
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGStraight::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	assert(length>0.);

	const MGVector& dir=direction();
	double tlen=length/dir.len();
	if(start){
		if(infinite_below())
			return;
		m_sparam-=tlen;
	}else{
		if(infinite_above())
			return;
		m_endparam+=tlen;
	}
	if(m_knotV){
		delete m_knotV;
		m_knotV=0;
	}
}

//Test if input parameter value is inside parameter range of the line.
bool MGStraight::in_range(double t)const{
	double error=MGTolerance::rc_zero()*m_direction.len();
				//error is maximum distance of 
				// two points that are assumed to be the same.
	double esave=MGTolerance::set_wc_zero(error);
	bool in=(m_sparam<=t && t<=m_endparam);
	MGTolerance::set_wc_zero(esave);
	return in;
}

//Test if this cure is co-planar with the 2nd curve curve2.
//MGPlane expression will be out to plane if this is co-planar.
//Function's return value is true if co-planar.
bool MGStraight::is_coplanar(const MGCurve& curve2, MGPlane& plane)const{
	const MGStraight* sl2=dynamic_cast<const MGStraight*>(&curve2);
	if(sl2){
		//When both are straight.
		if(direction().parallel(sl2->direction())){
			plane=MGPlane(*sl2,root_point());
			return true;
		}else{
			MGUnit_vector normal=direction()*sl2->direction();
			plane=MGPlane(normal,root_point());
			return plane.on(sl2->root_point());
		}
	}

	MGStraight line2;
	MGPosition point;
	int plkind=-1;
	const MGLBRep* lb=dynamic_cast<const MGLBRep*>(&curve2);
	if(lb)
		plkind=lb->planar(plane,line2,point);
	else{
		const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(&curve2);
		if(rlb)
			plkind=lb->planar(plane,line2,point);
	}
	if(plkind==0) return false;
	if(plkind==1){
		plane=MGPlane(*this,point);
		return true;
	}else if(plkind==2){
		if(!direction().parallel(line2.direction())) return false;
		plane=MGPlane(line2,root_point());
		return true;
	}else if(plkind==3){
		return plane.on(*this);
	}

	//When curve2 is neither MGLBRep nor MGRLBRep.
	if(!curve2.is_planar(plane)) return false;
	return plane.on(*this);
}

//Test if the input parameter t is the start point parameter or not.
bool MGStraight::is_startpoint_parameter(double t)const{
	if(m_sparam.minus_infinite())
		return false;

	if(m_endparam.plus_infinite()){
		return m_sparam==t;
	}

	return MGREqual_base(param_s(),t,param_span());
}

//Test if the input parameter t is the start point parameter or not.
bool MGStraight::is_endpoint_parameter(double t)const{
	if(m_endparam.plus_infinite())
		return false;

	if(m_sparam.minus_infinite()){
		return m_endparam==t;
	}

	return MGREqual_base(param_s(),t,param_span());
}

//Test if this cure is linear or not, that is, is straight or not.
//MGStraight expression will be out to straight if this is linear or not.
//Function's return value is true if linear.
bool MGStraight::is_linear(MGStraight& straight)const{
	straight=(*this);
	return true;
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGStraight::is_planar(MGPlane& plane)const{
	const MGVector& sldir=direction();
	double xx=sldir%mgX_UVEC; xx*=xx;
	double yy=sldir%mgY_UVEC; yy*=yy;
	double zz=sldir%mgZ_UVEC; zz*=zz;
	MGVector dir2;
	if(xx<=yy){
		if(xx<=zz)
			dir2=mgX_UVEC;
		else
			dir2=mgZ_UVEC;
	}else{
		if(yy<=zz)
			dir2=mgY_UVEC;
		else
			dir2=mgZ_UVEC;
	}
	plane=MGPlane(sldir,dir2,root_point());
	return true;
}

////////////isect with a curve.

// Straight と Curve の交点を求める。
MGCCisect_list MGStraight::isect(const MGCurve& curve)const{
	MGCCisect_list list=curve.isect(*this);
	return list.replace12();
}

// Straight と Straight の交点を求める。
MGCCisect_list MGStraight::isect(const MGStraight& st) const{
	MGCCisect p; MGCCisect_list list(this, &st);
	MGPSRELATION rel=relation(st,p);
	if(rel== MGPSREL_ISECT || rel==MGPSREL_COIN )
		list.append(p);
	return list;
}

//Compute intersections with MGRLBRep curve2.
MGCCisect_list MGStraight::isect(const MGRLBRep& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection with a ellipse.
MGCCisect_list MGStraight::isect(const MGEllipse& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection with a MGSurfCurve.
MGCCisect_list MGStraight::isect(const MGSurfCurve& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection with a MGBSumCurve.
MGCCisect_list MGStraight::isect(const MGBSumCurve& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list MGStraight::isect_with_noCompoSC(const MGSurfCurve& curve2)const{
	MGCCisect_list list=curve2.isect_noCompo(*this);
	return list.replace12();
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list MGStraight::isect_withC1LB(const MGLBRep& curve2)const{
	MGCCisect_list list=curve2.C1isect(*this);
	return list.replace12();
}

////////////isect with a surface.

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGSurface& surf) const{
	return surf.isect(*this);
}

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGPlane& surf) const{
	return surf.isect(*this);
}

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGSphere& surf) const{
	return surf.isect(*this);
}

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGCylinder& surf) const{
	return surf.isect(*this);
}

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGSBRep& surf) const{
	return surf.isect(*this);
}

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGRSBRep& surf) const{
	return surf.isect(*this);
}

//Intersection of Straight and aplane.
MGCSisect_list MGStraight::isect(const MGBSumSurf& surf) const{
	return surf.isect(*this);
}

MGCSisect_list MGStraight::isect(const MGFace& f)const{
	MGCSisect_list list;
	const MGSurface* srf=f.surface();
	if(!srf)
		return list;

	const MGBox& sbx=f.box();
	if(!sbx.crossing(*this))
		return list;

	list=srf->isectSl(*this,f.box_param());
	MGCSisect_list::CSiterator	i=list.begin(), iend=list.end(), i1;
	while(i!=iend){
		i1=i; i1++;
		if(!f.in_range((*i).param_surface()))
			list.removeAt(i);
		i=i1;
	}
	return list;
}

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
MGCParam_list MGStraight::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
)const{
	MGCParam_list list(this);

	double u=direction().ref(coordinate);
	if(MGMZero(u)){
		if(MGAEqual(f,m_root_point(coordinate)))
			list.append(m_sparam.value());
	}else{
		double p=f-root_point().ref(coordinate);
		double t=p/u;
		if(in_range(t)) list.append(t);
	}
	return list;
}

//isect2D returns parameter values of this(t) and l2(s)
// of the intersection point of both 2D straight lines.
// This and l2 are treated as infinite lines.
//Function's return value is:
// true if two lines are parallel(or one of the directin is zero)
// false if intersection was obtained.
bool MGStraight::isect2D(const MGStraight& l2,double& t,double& s)const{
	double x1, a1, y1, b1, x2, a2, y2, b2;
	x1=root_point().ref(0); a1=m_direction.ref(0);
	y1=root_point().ref(1); b1=m_direction.ref(1);
	x2=l2.root_point().ref(0); a2=l2.m_direction.ref(0);
	y2=l2.root_point().ref(1); b2=l2.m_direction.ref(1);
	double det = a1*b2-b1*a2;
	if(MGMZero(det))
		return true;

	double y1my2=y1-y2; double x2mx1=x2-x1;
	t=(a2*y1my2+b2*x2mx1)/det;
	s=(a1*y1my2+b1*x2mx1)/det;
	return false;
}

//Access to i-th element of knot.
//i=0, 1 and returns start or end parameter value of the ellipse.
double MGStraight::knot(size_t i) const{
	assert(i<=3);
	if(i<=1)
		return param_s();
	else
		return param_e();
}

//Returns the knot vector of the curve.
const MGKnotVector& MGStraight::knot_vector() const{
	if(!m_knotV){
		m_knotV=new MGKnotVector(2,2);
		(*m_knotV)(0)=(*m_knotV)(1)=param_s();
		(*m_knotV)(2)=(*m_knotV)(3)=param_e();
	}
	return *m_knotV;
}
MGKnotVector& MGStraight::knot_vector(){
	if(!m_knotV){
		m_knotV=new MGKnotVector(2,2);
		(*m_knotV)(0)=(*m_knotV)(1)=param_s();
		(*m_knotV)(2)=(*m_knotV)(3)=param_e();
	}
	return *m_knotV;
}
// パラメータ値が昇順であたえられた場合正の値で、降順の場合
// 負の値でパラメータ間の直線の代数的距離を返却する。
double MGStraight::length (double d1, double d2) const {
	return (range(d2) - range(d1))*m_direction.len();
}

// 自身の直線が有界の場合、その直線の距離を返却する。
// 非有界のとき－２を返却する。
double MGStraight::length () const {
	double d = -2.0;
	if (straight_type() == MGSTRAIGHT_SEGMENT)
		d=(m_endparam.value()-m_sparam.value())*m_direction.len();
	else if (straight_type() == MGSTRAIGHT_EMPTY) d = 0.0;
	return d;
}

// 直線上の与えられたパラメータで示される点から指定距離はなれた点
// のパラメータ値をかえす。
double MGStraight::length_param (double t,double len) const {
	return range(t+len/m_direction.len());
}

// 自身の直線に指定されたｌｉｍｉｔを付与する。
MGStraight& MGStraight::limit(const MGInterval& l) {
	if(m_knotV){delete m_knotV; m_knotV=0;}
	MGInterval p_limit = param_range(); // 直線範囲のパラメータを作成
	p_limit &= l;						// 指定Intervalとの積をとる
	m_sparam=p_limit.low();
	m_endparam=p_limit.high();
	update_mark();
	return *this;
}

//Compute sub straight line limitted by an box.
//This box's coordinates consist of world coordinates.
MGStraight& MGStraight::limit(const MGBox& box){
	MGInterval prange=param_range();
	MGCParam_list list;
	double bound, t;
	size_t dim=sdim();
	for(size_t i=0; i<dim; i++){
		bound=box.ref(i).low_point();
		list=isect_1D(bound, i);
		if(list.entries()){
			t=list.first();
			if(prange>>t){
				if(m_direction.ref(i)>=0.) prange.set_low_point(t);
				else                       prange.set_high_point(t);
			}
		}
		bound=box.ref(i).high_point();
		list=isect_1D(bound, i);
		if(list.entries()){
			t=list.first();
			if(prange>>t){
				if(m_direction.ref(i)>=0.) prange.set_high_point(t);
				else                       prange.set_low_point(t);
			}
		}
	}
	limit(prange);
	return *this;
}

//Compute nearest point on the line to the origin.
MGPosition MGStraight::nearest_to_origin() const{	
	double lensqr=m_direction.len(); lensqr*=lensqr;
	return MGPosition
		(root_point()-m_direction*((root_point()%m_direction)/lensqr));
}

// 直線の方向を反転する(方向ベクトルを逆向きにする）。
// 始終点があるときは入れ換える。
void MGStraight::negate(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	MGEReal save=m_sparam;
	m_sparam=-m_endparam; m_endparam=-save;
	m_direction = -m_direction;
}

//Obtain parameter value if this curve is negated by "negate()".
double MGStraight::negate_param(double t)const{
	return -t;
}

// 点が直線上にあるかを試験する。直線上にあれば，その点のパラメーター値を，
//　直線上になくても最近傍点のパラメーター値を返す。
bool MGStraight::on(
	const MGPosition& p,	 // 指定点       
	double& d				 // パラメータ値 
)const{
	bool on; double t;
	if(perp_point(p, t)) on=MGAZero(MGVector(p,eval(t)).len());
	else                 on=0;
	d=range(t);
	return on;
}

// 直線が平面上にあるか調べる。（平面上ならばtrue）
bool MGStraight::on(
	const MGPlane& pl	// Plane
)const{
	MGPosition uv(2);
	return direction().orthogonal(pl.normal()) && pl.on(root_point(), uv);
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGStraight::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
)const{
	MGStraight* sl=new MGStraight(*this);
	MGVector G(3,g);
	MGVector one(1);
	one(0)=G%m_direction; sl->m_direction=one;
	one(0)=G%m_root_point-g[3]; sl->m_root_point=one;
	return std::auto_ptr<MGCurve>(sl);
}

// 直線上の与えられたポイントにおけるパラメータ値を返す。
// If input point is not on the curve, return the nearest point on the
// curve.
double MGStraight::param(const MGPosition& p)const{
	// 垂線の足を求める。
	double d2;
	perp_point(p,d2);

	// パラメータ値を返す。
	return range(d2);
}

//Obtain parameter space error.
double MGStraight::param_error()const{
	MGEReal len=m_endparam-m_sparam;
	if(len.finite())
		return len.value()*MGTolerance::rc_zero();
	else{
		double vlen=m_direction.len();
		return vlen*MGTolerance::rc_zero();
	}
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance. For straight, the knots are start and end points.
double MGStraight::param_normalize(double t) const{
	double plen;
	if(straight_type()==MGSTRAIGHT_SEGMENT) plen=param_span();
	else plen=1.;
	double tnew;
	if(straight_type()==MGSTRAIGHT_SEGMENT
		|| straight_type()==MGSTRAIGHT_HALF_LIMIT){
		tnew=param_s();
		if(MGRZero2(t-tnew,plen)) return tnew;
	}
	if(straight_type()==MGSTRAIGHT_SEGMENT){
		tnew=param_e();
		if(MGRZero2(t-tnew,plen)) return tnew;
	}
	return t;
}

//Return parameter range of the curve(パラメータ範囲を返す)
MGInterval MGStraight::param_range()const{
	return MGInterval(m_sparam, m_endparam);
}

//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGStraight* MGStraight::part(double t1, double t2, int multiple)const{
	assert(t2>=t1);
	MGStraight* sl=new MGStraight(*this);
	sl->limit(MGInterval(t1,t2));
	return sl;
}

// 与点から直線への垂線の足のパラメータ値を返却する。
//Return the foot of the  straight line that is perpendicular to this line.
//Function's return value is parameter value of this straight line,
//may ***NOT*** be in_range.
double MGStraight::perp_param(
	const MGPosition& point		// 与点
)const{
	MGVector v2(point,root_point());// 始点から与点へのベクトル			
	double lensqr=m_direction.len(); lensqr*=lensqr;
	return v2%m_direction/lensqr;	// 与点から直線への垂線の足と始点との距離
}

// 与えられたポイントから曲線への垂線の足、そのパラメータ値を返却する。
// Function's return value is if point is obtained(1) or not(0)
int MGStraight::perp_point(
	const MGPosition& point,	// 指定点
	double& d1,					// 垂線の足のパラメータ値
	const double* d2			// guess parameter value of d1.
)const{
	d1=perp_param(point);	
	return in_range(d1);
}
	
// 与ポイントから直線へ下ろした垂線の足の，直線のパラメータ値を
// すべて求める。
MGCParam_list MGStraight::perps(
	const MGPosition& point	// 与ポイント
)const{
	MGCParam_list tlist(this);
	double l =perp_param(point);	
	if(in_range(l)) tlist.append(l);
	return tlist;
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv2's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list MGStraight::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGStraight::perps(
	const MGStraight& l2		//The second curve
)const{
	if(MGMZero(m_direction.sangle(l2.m_direction)))// 平行
		return relation_parallel(l2);
	MGPosition P;
	perp_guess(1.,0.,l2,1.,0.,0.,0.,P);
	MGPosition_list list;
	if(in_range(P.ref(0)) && l2.in_range(P.ref(1)))
		list.append(P);
	return list;
}
MGPosition_list MGStraight::perps(
	const MGRLBRep& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGStraight::perps(
	const MGEllipse& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGStraight::perps(
	const MGSurfCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGStraight::perps(
	const MGBSumCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

// 入力パラメータをパラメータ範囲で丸めて返却する。
double MGStraight::range(double d) const {
	MGInterval limit = param_range();
	if (straight_type() == MGSTRAIGHT_EMPTY) d=0.0;
	else if (limit > d)						 d=m_sparam.value();
	else if (limit < d)						 d=m_endparam.value();
	return d;
}

// Compute two straight lines relationship.
// (自身と与えられた直線の関係を調べる。)
MGPSRELATION MGStraight::relation_parallel (
	const MGStraight& s,
	MGCCisect& ip ) const{
	double t1=0.,t2=0.; MGPSRELATION rl;
	if(m_direction.parallel(s.root_point()-root_point())){
		rl = MGPSREL_COIN;
		if(straight_type()==MGSTRAIGHT_UNLIMIT){
			t2=s.param_s(); t1=param(s.root_point());
		}else if(s.straight_type()==MGSTRAIGHT_UNLIMIT){
			t1=param_s(); t2=s.param(root_point());
		}else if(straight_type()==MGSTRAIGHT_HALF_LIMIT &&
			s.straight_type()==MGSTRAIGHT_HALF_LIMIT){
			MGPosition sp1,sp2;
			if(m_sparam.finite()) sp1=start_point();
			else sp1=end_point();
			if(s.m_sparam.finite()) sp2=s.start_point();
			else sp2=s.end_point();
			if(s.on(sp1,t2)) t1=param_s();
			else if(on(sp2,t1)) t2=s.param_s();
			else {rl=MGPSREL_PARALLEL; t1=t2=0.;}
		}else if(straight_type()==MGSTRAIGHT_SEGMENT){
			if(s.on(start_point(),t2)) t1=param_s();	// 一致
			else if(s.on(end_point(),t2)) t1=param_e();	// 一致 
			else {rl = MGPSREL_PARALLEL; t2=0.;}		// 平行
		}else{		// s is MGSTRAIGHT_SEGMENT.
			if(on(s.start_point(),t1)) t2=s.param_s();	// 一致
			else if(on(s.end_point(),t1)) t2=s.param_e();	// 一致 
			else {rl = MGPSREL_PARALLEL; t1=0.;}		// 平行
		}
		ip=MGCCisect(eval(t1),t1,t2,MGCCREL_COIN);
	}
	else rl = MGPSREL_PARALLEL;
	return rl;
}

// 自身と与えられた直線の関係を調べる。
//Return two straight line's relationship.
//MGPSREL_UNKNOWN,     //Unknown. 不明
//MGPSREL_TORSION,     //Torsion. ねじれ
//MGPSREL_ISECT,       //Intersection. 交差
//MGPSREL_PARALLEL,    //Parallel. 平行
//MGPSREL_COIN,        //Coincidence. 包含
//MGPSREL_VIRTUAL_ISECT//Virtually intersection
//                     (Intersection at extended line).
//                     指定直線の延長上の点での交差
MGPSRELATION MGStraight::relation (
	const MGStraight& s,
	MGCCisect& ip
)const{	
	double t1=0.,t2=0.; MGPSRELATION rl;
	ip=MGCCisect();
	if(straight_type()==MGSTRAIGHT_EMPTY ||
		s.straight_type()==MGSTRAIGHT_EMPTY) return MGPSREL_TORSION;

	// ２直線の関係のチェック、交点なしかどうか。
	if(s.sdim()==2 && sdim()==2){
		if(isect2D(s, t1, t2)) return relation_parallel(s,ip);
	}else if(MGMZero(m_direction.sangle(s.m_direction))){// 平行
		return relation_parallel(s,ip);
	}else{
		// 二直線から自身の直線を含む平面を作成する。
		MGPlane pl(m_direction, s.m_direction, root_point());
		//pl's normal is normal to both m_direction and s.m_direction.
		MGMatrix mat(3); mat.set_axis(pl.normal(),2);
		MGStraight sl1=(*this)*mat; //sl1 is 2D straightline on x-y plane.
		MGStraight sl2=s*mat;		//sl2 is 2D straightline on x-y plane.
		if(sl1.isect2D(sl2,t1,t2)){	// パラメータ値の計算
			return MGPSREL_TORSION;
		}
	}

	MGPosition point(root_point() + t1*m_direction); //Intersection Point.
	MGCCRELATION ccrel; 
	// 交点が線分上か延長線上にあるかを調べる。
	if(in_range(t1) && s.in_range(t2)) {
		ccrel=MGCCREL_ISECT;
		rl=MGPSREL_ISECT;
	}else{
		ccrel=MGCCREL_UNKNOWN;
		rl=MGPSREL_VIRTUAL_ISECT;
	}
	ip=MGCCisect(point,t1,t2,ccrel);

	return rl; //関係をかえす。
}

// 自身と与えられた平面の関係を調べる。
MGPSRELATION MGStraight::relation (
		const MGPlane & pl,
		MGCSisect & ip ) const{
	MGPSRELATION rl; MGCSRELATION csrel;
	// 直線の方向ベクトルと平面の法線ベクトルの内積を求める。
	double cross = pl.normal()%m_direction;
	if(MGMZero(cross)){	// 直線と平面が平行、又は直線が平面上
		if (MGAZero(pl.distance(root_point()))){	// 直線が平面上
			rl = MGPSREL_COIN; csrel=MGCSREL_COIN;
		}else{			// 直線と平面が平行
			rl = MGPSREL_PARALLEL; csrel=MGCSREL_UNKNOWN;
		}
		ip=MGCSisect(root_point(), 0., MGPosition(0.,0.),csrel);
	}else{				// 交差するとき
		// 直線上のパラメータ値を求める。
		double t=(pl.distance()-pl.normal()%root_point())/cross;
		// 交点を求める。
		MGPosition point = MGPosition(root_point() + t*m_direction);
		MGPosition uv(2);
		pl.on(point,uv);	//Compute plane's parameter value of point.
		ip = MGCSisect(point,t,uv,MGCSREL_IN);

		// 交点が直線上にあるかしらべる。
		if(in_range(t))	 // 直線上にあるとき
			rl = MGPSREL_ISECT;
		else			 // 直線上にないとき
			rl = MGPSREL_VIRTUAL_ISECT;
	}
	return rl;
}

//Compute parallel range of two straight lines.
//Two straight line this and l2 must be parallel.
MGPosition_list MGStraight::relation_parallel(const MGStraight& l2) const
//Function's return value MGPosition_list list is:
// list.entries() is 0(when no paralle part) or 2(when there is a parallel
// part). When list.entries() is 2, let their entries be P1 and P2.
// Then from P1(0) to P2(0) is the range of this straight line.
// From P1(1) to P2(1) is the range of line2 straight line.
{
	assert(m_direction.parallel(l2.m_direction));

	MGPosition_list list;
	double s1=perp_param(l2.start_point()), s2=perp_param(l2.end_point());
	double t1=l2.perp_param(start_point()), t2=l2.perp_param(end_point());
	//s1, s2: perp point on this from l2 start and end point
	//t1, t2: perp point on l2 from this start and end point
	MGPosition Ps(param_s(),t1), Pe(param_e(),t2);
	MGPosition Qs(s1,l2.param_s()), Qe(s2,l2.param_e());

	if(straight_type()==MGSTRAIGHT_EMPTY ||
		l2.straight_type()==MGSTRAIGHT_EMPTY) return list;

	if(straight_type()==MGSTRAIGHT_SEGMENT){		//When this is segment.
		if(l2.in_range(t1)){
			list.append(Ps);
			if(l2.in_range(t2)) list.append(Pe);
			else{
				if(in_range(s1)) list.append(Qs);
				else list.append(Qe);
			}
		}else if(l2.in_range(t2)){
			if(in_range(s1)) list.append(Qs);
			else list.append(Qe);
			list.append(Pe);
		}else if(t1*t2<0.){
			list.append(Qs); list.append(Qe);
		}
	}
	else if(straight_type()==MGSTRAIGHT_UNLIMIT){	//When this is unlimit.
		list.append(Qs); list.append(Qe);
	}else if(l2.straight_type()==MGSTRAIGHT_UNLIMIT){//When l2 is unlimt.
		list.append(Ps); list.append(Pe);
	}else{										//When this is half limit.
		if(l2.straight_type()==MGSTRAIGHT_SEGMENT){	//l2 is segment.
			if(in_range(s1)){
				list.append(Qs);
				if(in_range(s2)) list.append(Qe);
				else list.append(Ps);
			}else if(in_range(s2)){
				list.append(Ps); list.append(Qe);
			}
		}else{									//When both are half_limit.
			if(in_range(s1)){
				if(l2.in_range(t1)){
					list.append(Ps); list.append(Qs);
				}else{
					list.append(Qs); list.append(Qe);
				}
			}else if(l2.in_range(t1)){
				list.append(Ps); list.append(Pe);
			}
		}
	}
	return list;
}

//Function to avoid m_direction.len()=zero.
//The line's m_direction is set as a unit vector and m_endparam
//is set to zero.
void MGStraight::save_length_zero(){
	double dlen=m_direction.len();
	if(MGMZero(dlen)){
		m_direction=m_direction.normalize();
		m_endparam=m_sparam+dlen;
		if(m_knotV){delete m_knotV; m_knotV=0;}
	}
}

//Return space dimension
size_t MGStraight::sdim() const{	
	size_t dim1=m_root_point.sdim();
	if(dim1==0) return 0;
	size_t dim2=m_direction.sdim();
	if(dim1<dim2) dim1=dim2;
	return dim1;
}

// 直線のタイプ，方向ベクトル，始点を指定して直線を生成する。
//Straight from straight line type, direction vector, and an origin.
//Construct a straight and replce this with it.
//This fuction does not convert input vec to a unit vector.
//If you like the conversion, use MGStraight() constructor.
MGStraight& MGStraight::set_straight(
	MGSTRAIGHT_TYPE type,				//Type
	const MGVector& vec,				//Direction
	const MGPosition& Q					//Origin
){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	m_root_point=Q;
	m_direction=vec;
	m_sparam=0.;
	m_endparam=MGEReal(MGINFINITE_PLUS);

	if (type==MGSTRAIGHT_SEGMENT) m_endparam=1.;
	else if(type==MGSTRAIGHT_EMPTY) m_endparam=-1.;
	else if(type==MGSTRAIGHT_UNLIMIT) m_sparam=MGEReal(MGINFINITE_MINUS);
	size_t dimv=vec.sdim(), dimp=Q.sdim();
	if(dimv<dimp) m_direction=MGVector(dimp,vec);
	else if(dimv>dimp) m_root_point=MGPosition(dimv,Q);
	save_length_zero();
	update_mark();
	return *this;
}

// 自身の直線の始点を返却する。
//Return start(root) point of the straight.
MGPosition MGStraight::start_point() const{
	return eval_position(m_sparam.value());
}

// 直線のタイプを返却する。
//Return the straight line type.
MGSTRAIGHT_TYPE MGStraight::straight_type() const{
	if(m_sparam>m_endparam) return MGSTRAIGHT_EMPTY;
	else if(m_sparam.finite()){
		if(m_endparam.finite()) return MGSTRAIGHT_SEGMENT;
		else return MGSTRAIGHT_HALF_LIMIT;
	}else{
		if(m_endparam.finite()) return MGSTRAIGHT_HALF_LIMIT;
		else return MGSTRAIGHT_UNLIMIT;
	}
}

// 自身の直線からｌｉｍｉｔを取り除く。
MGCurve& MGStraight::unlimit(){
	m_sparam=MGEReal(MGINFINITE_MINUS);
	m_endparam=MGEReal(MGINFINITE_PLUS);
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGStraight::unlimit_end(){
	m_endparam=MGEReal(MGINFINITE_PLUS);
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGStraight::unlimit_start(){
	m_sparam=MGEReal(MGINFINITE_MINUS);
	update_mark();
	return *this;
}

//
// 演算子の多重定義
//

//Assignment.
//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGStraight& MGStraight::operator=(const MGStraight& sl2){
	if(this==&sl2)
		return *this;

	MGCurve::operator=(sl2);
	m_direction=sl2.m_direction;
	m_endparam=sl2.m_endparam;
	if(m_knotV){
		delete m_knotV; m_knotV=0;
	}
	m_root_point=sl2.m_root_point;
	m_sparam=sl2.m_sparam;
	return *this;
}
MGStraight& MGStraight::operator=(const MGGel& gel2){
	const MGStraight* gel2_is_this=dynamic_cast<const MGStraight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

// 直線の平行移動を行いオブジェクトを生成する。
MGStraight MGStraight::operator+(const MGVector& vec) const{
	MGStraight s = *this;
	s += vec;
	return s;
}
MGStraight operator+(const MGVector& v, const MGStraight& sl){
	return sl+v;
}

// 直線の平行移動を行い自身の直線とする。
MGStraight& MGStraight::operator+=(const MGVector & v){
	m_root_point += v;
	if(m_box) (*m_box)+=v;
	return *this;
}

// 直線の逆方向に平行移動を行いオブジェクトを生成する。
MGStraight MGStraight::operator-(const MGVector & v)const{
	MGStraight s = *this;
	s -= v;
	return s;
}

//直線の逆方向に平行移動を行い自身の直線とする。
MGStraight& MGStraight::operator-= (const MGVector & v){
	m_root_point -= v;
	if(m_box) (*m_box)-=v;
	return *this;
}

// 与えられたスケールで直線の変換を行いオブジェクトを生成する。
//Generate scaled straight.
MGStraight MGStraight::operator*(double scale) const{
	MGStraight sl(*this);
	sl *=scale;
	return sl;
}

// 与えられたスケールで直線の変換を行いオブジェクトを生成する。
//Generate scaled straight.
MGStraight operator*(double scale, const MGStraight& sl){
	return sl*scale;
}

// 与えられたスケールで直線の変換を行い自身の直線とする。
//Scaling the straight.
MGStraight& MGStraight::operator*=(double scale){
	m_root_point *=scale;
	m_direction *=scale;
	save_length_zero();
	update_mark();
	return *this;
}

// 与えられた変換で直線の変換を行いオブジェクトを生成する。
MGStraight MGStraight::operator*(const MGMatrix& t)const{
	MGStraight s = *this;
	s *= t;
	return s;
}

// 与えられた変換で直線の変換を行い自身の直線とする。
MGStraight& MGStraight::operator*=( const MGMatrix& t ){
	m_root_point *= t;
	m_direction = m_direction*t;
	save_length_zero();
	update_mark();
	return *this;
}

// 与えられた変換で直線のトランスフォームを行いオブジェクトを生成する。
MGStraight MGStraight::operator*(const MGTransf& t)const{
	MGStraight s = *this;
	s *= t;
	return s;
}

// 与えられた変換で直線のトランスフォームを行い自身の直線とする。
MGStraight& MGStraight::operator*=(const MGTransf& t){
	m_root_point *= t;
	m_direction *=t.affine();
	save_length_zero();
	update_mark();
	return *this;
}

//
// 論理演算子の多重定義
bool MGStraight::operator==(const MGStraight& sl2)const{
	if(sdim()==0 && sl2.sdim()==0)
		return 1;
	 return (m_root_point == sl2.m_root_point 
		&& m_direction == sl2.m_direction
		&& m_sparam==sl2.m_sparam
		&& m_endparam==sl2.m_endparam);
}

bool MGStraight::operator<(const MGStraight& gel2)const{
	return m_root_point.len()<gel2.m_root_point.len();
}
bool MGStraight::operator==(const MGGel& gel2)const{
	const MGStraight* gel2_is_this=dynamic_cast<const MGStraight*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGStraight::operator<(const MGGel& gel2)const{
	const MGStraight* gel2_is_this=dynamic_cast<const MGStraight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}
