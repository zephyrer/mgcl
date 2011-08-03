/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/BSumCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/SurfCurve.h"
#include "mg/nlbit.h"
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

//Copy constructor
MGSurfCurve::MGSurfCurve(const MGSurfCurve& sc)
:MGCurve(sc),m_surface(sc.m_surface), m_curve(sc.m_curve){;}

//Default constructor
MGSurfCurve::MGSurfCurve(const MGSurface& srf, const MGCurve& crv)
:MGCurve(crv),m_surface(&srf),m_curve(crv, crv.param_range()){
}

//Default constructor
MGSurfCurve::MGSurfCurve(
	const MGFSurface& srf, const MGCurve& crv
):MGCurve(crv),m_surface(srf.get_surface_pointer()),
m_curve(crv, crv.param_range()){
}

//Test if m_curve is MGCompositeCurve. If composite, return
//the pointer. If not, return null.
const MGCompositeCurve* MGSurfCurve::base_composite()const{
	const MGCompositeCurve* compo
		=dynamic_cast<const MGCompositeCurve*>(m_curve.base_curve());
	return compo;
}

//Changing this object's space dimension.
MGSurfCurve& MGSurfCurve::change_dimension(
	size_t dim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2		// Source order of this object.
){ 
	assert(false);//This cannot be used.
	return *this;
}

//Return minimum box that includes whole of the curve.
//曲線部分を囲むボックスを返す。
MGBox* MGSurfCurve::compute_box()const{
	MGBox cbox=box_limitted(param_range());
	return new MGBox(cbox);
}

//Return minimum box that includes the curve of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox MGSurfCurve::box_limitted(
	const MGInterval& l // Parameter Range of the curve.
)const{
	//uv-curveを囲むボックス(uv-interval)を内包する
	//m_surfaceのボックス。
	return m_surface->box_limitted(m_curve.box_limitted(l));
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGSurfCurve* MGSurfCurve::clone()const{
	return new MGSurfCurve(*this);
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
// とりあえず3次微分まで展開して計算しています。
// 続きは後で考えましょう。
MGVector MGSurfCurve::eval(
	double t,				// Parameter value.
	size_t nderiv,		// Order of Derivative.
	int left			//Left continuous(left=true)
						//or right continuous(left=false).
)const{
//	assert(nderiv < 4);
	const size_t ncd = sdim();
	MGVector ret(ncd, 0.); //零ベクトルで初期化
	MGPosition prm_uv= m_curve.eval(t, 0, left);//母曲面のuvパラメータを求める。
	switch (nderiv){
	case 0:{
		ret = m_surface->eval(prm_uv, 0, 0);
		break;
		   }
	case 1:{
		MGVector dudt = m_curve.eval(t, 1, left);
		ret = m_surface->eval(prm_uv, 1, 0) * dudt.ref(0)
			+ m_surface->eval(prm_uv, 0, 1) * dudt.ref(1);
		break;
		   }

	case 2:	{
		MGVector dudt1 = m_curve.eval(t, 1, left);
		double dudt1_0 = dudt1.ref(0);
		double dudt1_1 = dudt1.ref(1);
		MGVector dudt2 = m_curve.eval(t, 2, left);
		ret = m_surface->eval(prm_uv, 2, 0) * dudt1_0*dudt1_0
			+ m_surface->eval(prm_uv, 1, 1) * dudt1_0*dudt1_1 *2.
			+ m_surface->eval(prm_uv, 1, 0) * dudt2.ref(0)
			+ m_surface->eval(prm_uv, 1, 0) * dudt2.ref(1)
			+ m_surface->eval(prm_uv, 0, 2) * dudt1_1*dudt1_1;
		break;
			}
	case 3:{
		MGVector dudt1  = m_curve.eval(t, 1, left);
		MGVector dudt2 = m_curve.eval(t, 2, left);
		MGVector dudt3 = m_curve.eval(t, 3, left);
		double dudt1_0 = dudt1.ref(0);
		double dudt1_1 = dudt1.ref(1);
		double dudt2_0 = dudt2.ref(0);
		double dudt2_1 = dudt2.ref(1);

		ret = m_surface->eval(prm_uv, 3, 0) * dudt1_0*dudt1_0*dudt1_0
			+ m_surface->eval(prm_uv, 2, 1) * dudt1_0*dudt1_0*dudt1_1 *3.
			+ m_surface->eval(prm_uv, 1, 2) * dudt1_0*dudt1_1*dudt1_1 *3.
			+ m_surface->eval(prm_uv, 0, 3) * dudt1_1*dudt1_1*dudt1_1
			+ m_surface->eval(prm_uv, 2, 0) * dudt2_0*dudt1_0 * 3.
			+ m_surface->eval(prm_uv, 1, 1) * dudt2_0*dudt1_1 * 3.
			+ m_surface->eval(prm_uv, 1, 1) * dudt2_1*dudt1_0 * 3.
			+ m_surface->eval(prm_uv, 0, 2) * dudt2_1*dudt1_1 * 3.
			+ m_surface->eval(prm_uv, 1, 0) * dudt3.ref(0)
			+ m_surface->eval(prm_uv, 1, 0) * dudt3.ref(1);
		break;
		   }
	default:
		{ //4次微分以上の場合は零ベクトルを返す。(初期値のまま)
		break;
		}
	}

	return ret;
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGSurfCurve::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	double ts=param_s(), te=param_e();
	double t=start ? ts:te;
	MGVector V1=eval(t,1);
	MGVector V2=m_curve.eval(t,1);
	double tan_ratio=V1.len()/V2.len();
	length/=tan_ratio;
	if(start)
		length*=-1.;

	m_curve.extend(length,start);
	double prm[4]={m_surface->param_s_u(),m_surface->param_s_v(),
		m_surface->param_e_u(),m_surface->param_e_v()};
	for(size_t i=0; i<4; i++){
		MGCParam_list plist=m_curve.isect_1D(prm[i],i%2);
		if(plist.size()){
			if(start)
				m_curve.limit(MGInterval(plist.back(),te));
			else
				m_curve.limit(MGInterval(ts,plist.front()));
		}
	}
}

//Test if input parameter value is inside parameter range of the line.
bool MGSurfCurve::in_range(double t) const{
	return m_curve.in_range(t);
}

//Provide divide number of curve span for function intersect.
size_t MGSurfCurve::intersect_dnum() const{
	size_t n1=m_curve.intersect_dnum(), n2,n3;
	const MGRSBRep* surf_is_MGRSBRep=dynamic_cast<const MGRSBRep*>(m_surface);

	const MGBox& uvbox=m_curve.box();
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	const MGKnotVector& tu=m_surface->knot_vector_u();
	int ordr=tu.order();
	size_t i1=tu.locate(u0), i2=tu.locate(u1,1);
	if(i2<i1)
		i2=i1;
	if(surf_is_MGRSBRep)
		n2=(i2-i1+1)*ordr;
	else
		n2=(i2-i1+1)*(ordr-1);

	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	const MGKnotVector& tv=m_surface->knot_vector_v();
	ordr=tv.order();
	i1=tv.locate(v0); i2=tv.locate(v1,1);
	if(i2<i1)
		i2=i1;
	if(surf_is_MGRSBRep)
		n3=(i2-i1+1)*ordr;
	else
		n3=(i2-i1+1)*(ordr-1);

	//uvカーブのintersect_dnum() + Surfaceのintersect_dnum()のうち
	//最大のもので代用する。
	if(n1<n2){
		if(n2>n3)
			return n2;
		else
			return n3;
	}else if(n3>n1)
		return n3;
	else
		return n1;
}

//Intersection with a curve.
MGCCisect_list MGSurfCurve::isect(const MGCurve& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//Intersection with a MGStraight.
MGCCisect_list MGSurfCurve::isect(const MGStraight& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//Intersection with a MGStraight.
MGCCisect_list MGSurfCurve::isect(const MGRLBRep& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//Intersection with a MGStraight.
MGCCisect_list MGSurfCurve::isect(const MGEllipse& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//Intersection with a MGStraight.
MGCCisect_list MGSurfCurve::isect(const MGLBRep& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//Intersection with a MGStraight.
MGCCisect_list MGSurfCurve::isect(const MGSurfCurve& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//Intersection with a MGStraight.
MGCCisect_list MGSurfCurve::isect(const MGBSumCurve& curve2) const{
	MGCCisect_list list(this, &curve2);
	isect_of_each(curve2,list);
	return list;
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list MGSurfCurve::isect_noCompo(const MGCurve& curve2)const{
	MGCCisect_list list=curve2.isect_with_noCompoSC(*this);
	return list.replace12();
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list MGSurfCurve::isect_noCompo(const MGSurfCurve& curve2)const{
	MGCCisect_list list(this, &curve2);
	if(curve2.surface()==surface()){
		MGCCisect_list pisects=m_curve.isect(curve2.m_curve);
		size_t n=pisects.size();
		if(n==0)
			return list;
		MGCCisect_list::iterator i=pisects.begin(), ie=pisects.end();
		for(; i!=ie; i++){
			MGCCisect& isi=*i;
			double t1=isi.param1(), t2=isi.param2();
			MGPosition P=eval(t1);
			if(P!=curve2.eval(t2)){
				MGPosition t12(t1,t2);
				perp_guess(1.,-1.,curve2,1.,-1.,t1,t2,t12);
				t1=t12[0]; t2=t12[1];
				P=eval(t1);
			}
			list.append(MGCCisect(P,t1,t2));
		}
		return list;
	}

	const MGCompositeCurve* compo2=curve2.base_composite();
	if(compo2){
		double ts=curve2.m_curve.param_s(), te=curve2.m_curve.param_e();
		MGCompositeCurve::const_iterator
			icurve=compo2->begin(), endcurve=compo2->end();
		const MGSurface& srf2=*(curve2.m_surface);
		for(; icurve!=endcurve; endcurve++){
			const MGCurve& crvi=**icurve;
			if(crvi.param_e()<=ts || te<=crvi.param_s())
				continue;
			MGSurfCurve scrv2(srf2,crvi);
			MGCCisect_list listi=intersect(scrv2);
			MGCCisect_list::iterator j=listi.begin(), je=listi.end();
			for(; j!=je; j++){
				if(in_range(j->param1()))
					list.append(*j);
			}
		}
		return list;
	}

	return intersect(curve2);
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list MGSurfCurve::isect_noCompo(const MGStraight& sline) const{
	MGCCisect_list list(this, &sline);

	double tmid=(param_s()+param_e())*0.5;
	MGVector N=m_surface->normal(m_curve.eval(tmid));
	MGPlane plane(sline.direction(), N, sline.root_point());
	//The plane that includes sline.
	MGCSisect_list cslist=isect(plane);//直線を含む平面と面上線との交点を求める。

	MGPosition p; double t1,t2;
	MGCSisect_list::CSiterator i;
	for(i=cslist.begin(); i!=cslist.end(); i++){
		t1=(*i).param_curve(); //交点のSurfCurve上のパラメータ値を求める。
		p=(*i).point();//交点の座標値を求める(eval(t1)よりはこちらのほうが良い？)
		if(sline.on(p,t2)) //交点が直線sline上に乗っていれば面上線と直線の交点になる。
			list.append(MGCCisect(p,t1,t2));
	}
	return list;
}

//isect of each elements of this m_curve,
//if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void MGSurfCurve::isect_of_each(
	const MGCurve& curve2,	//The isect objective curve.
	MGCCisect_list& list	//Obtained isect will be appended.
)const{
	const MGCompositeCurve* compo=base_composite();
	if(!compo){
		list=isect_noCompo(curve2);
		return;
	}

	MGCompositeCurve::const_iterator
		icurve=compo->begin(), endcurve=compo->end();
	double ts=m_curve.param_s(), te=m_curve.param_e();
	for(; icurve!=endcurve; endcurve++){
		const MGCurve& crvi=**icurve;
		if(crvi.param_e()<=ts || te<=crvi.param_s())
			continue;
		MGSurfCurve scrv(*m_surface,crvi);
		MGCCisect_list listi=scrv.isect_noCompo(curve2);
		MGCCisect_list::iterator j=listi.begin(), je=listi.end();
		for(; j!=je; j++){
			if(in_range(j->param1()))
				list.append(*j);
		}
	}
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGSurface& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGPlane& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGSphere& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGCylinder& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGSBRep& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGRSBRep& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//Intersection with a Surface
MGCSisect_list MGSurfCurve::isect(const MGBSumSurf& surf)const{
	MGCSisect_list list(this,&surf);
	isect_of_each(surf,list);
	return list;
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCSisect_list MGSurfCurve::isect_noCompo(const MGSurface& surf)const{
	return surf.isect_with_noCompoSC(*this);
}

class MGSurfCurvePlaneDist{
	const MGSurfCurve* m_curve;
	const MGPlane* m_plane;
public:
	MGSurfCurvePlaneDist(const MGSurfCurve* scrv, const MGPlane* plane)
		:m_curve(scrv), m_plane(plane){;};
	double operator()(double t)const{
		return m_plane->distance(m_curve->eval(t));};
};

// 面上線と平面の交点を求める。
MGCSisect_list MGSurfCurve::isect_noCompo(const MGPlane& pl) const{
	size_t ndiv=intersect_dnum();
	double t0=param_s(), t1=param_e();
	double delta=(t1-t0)/double(ndiv);

	double error=MGTolerance::wc_zero();
	MGCSisect_list list(this,&pl);

	MGPosition Ppre, Paft;
	double tpre=t1,taft=t0;
	double dpre=pl.distance(Ppre=eval(tpre));
	double daft=pl.distance(Paft=eval(taft));
	if(fabs(daft)<error*.5) list.append(Paft,t0,pl.uv(Paft));
	if(fabs(dpre)<error*.5) list.append(Ppre,t1,pl.uv(Ppre));

	//Prepare for mgNlbit.
	MGSurfCurvePlaneDist pdist(this,&pl);

	//Iterate by checking singned distance dpre and daft from the plane.
	//When dpre and daft have different signs, an intersection point must 
	//lie between tpre and taft.
	size_t i=0;
	while(i<=ndiv){
		dpre=daft; tpre=taft; Ppre=Paft;
		while(fabs(dpre)<=error){
			list.append(Ppre,tpre,pl.uv(Ppre));
			if(i>=ndiv) break;
			i++; tpre=t0+delta*double(i);
			Ppre=eval(tpre); dpre=pl.distance(Ppre);
		}
		if(i>=ndiv) break;
		i++; taft=t0+delta*double(i);
		Paft=eval(taft); daft=pl.distance(Paft);
		if(fabs(daft)<=error) continue;
		else if(dpre*daft<0.){
		//Now there exists a solution between tpre and taft.
			int ier;
			double x=mgNlbit(pdist, tpre,taft, error, 20, ier);
			MGPosition Pt=eval(x);
			list.append(Pt,x,pl.uv(Pt));
		}
	}

	return list;
}

//isect of each elements of this m_curve,
//if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void MGSurfCurve::isect_of_each(
	const MGSurface& surf,	//The isect objective surface.
	MGCSisect_list& list	//Obtained isect will be appended.
)const{
	const MGCompositeCurve* compo=base_composite();
	if(!compo){
		list=surf.isect_with_noCompoSC(*this);
		return;
	}

	MGCompositeCurve::const_iterator
		icurve=compo->begin(), endcurve=compo->end();
	double ts=m_curve.param_s(), te=m_curve.param_e();
	for(; icurve!=endcurve; endcurve++){
		const MGCurve& crvi=**icurve;
		if(crvi.param_e()<=ts || te<=crvi.param_s())
			continue;
		MGSurfCurve scrv(*m_surface,crvi);
		MGCSisect_list listi=scrv.isect_noCompo(surf);
		MGCSisect_list::iterator j=listi.begin(), je=listi.end();
		for(; j!=je; j++){
			if(in_range(j->param_curve()))
				list.append(*j);
		}
	}
}

//Access to i-th element of knot.
double MGSurfCurve::knot(size_t i) const{ return m_curve.knot(i);}

// Return ending parameter value.
double MGSurfCurve::param_e() const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_e();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGSurfCurve::param_normalize(double t) const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_normalize(t);
}

//Return parameter range of the curve(パラメータ範囲を返す)
MGInterval MGSurfCurve::param_range() const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_range();
}

// Return starting parameter value.
double MGSurfCurve::param_s() const{
	// uvカーブのパラメータが面上線のパラメータになる。
	return m_curve.param_s();
}
	
//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGSurfCurve* MGSurfCurve::part(double t1, double t2, int multiple) const{
	assert(t2-t1>param_error());
	MGSurfCurve* ret = new MGSurfCurve(*this);
	if(multiple){
		MGCurve* crv=m_curve.part(t1,t2,multiple);
		ret->m_curve=MGTrimmedCurve(*crv,crv->param_range());
		delete crv;
	}else{
		ret->m_curve.limit(MGInterval(t1,t2));
	}
	ret->update_mark();
	return ret;
}

//Round t into curve's parameter range.
// 入力パラメータをパラメータ範囲でまるめて返却する。
double MGSurfCurve::range(double t) const{
	return m_curve.range(t);
}

//Return space dimension
size_t MGSurfCurve::sdim() const{
	return m_surface->sdim();
}

//Return curve type(曲線のタイプを返す)
MGCURVE_TYPE MGSurfCurve::type() const{
	return MGCURVE_SURFACE;
}

//Operator overload(演算子多重定義)

//Assignment.
//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGSurfCurve& MGSurfCurve::operator=(const MGSurfCurve& gel2){
	if(this==&gel2)
		return *this;

	MGCurve::operator=(gel2);
	m_curve=gel2.m_curve;
	m_surface=gel2.m_surface;
	return *this;
}
MGSurfCurve& MGSurfCurve::operator=(const MGGel& gel2){
	const MGSurfCurve* gel2_is_this=dynamic_cast<const MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Logical operator overload(論理演算子多重定義)
//Test if two curves are equal.
// 与曲線と自身が等しいかの比較判定を行う。
// 等しくなるのは自分自身と、自分自身から作ったコピーのみ。
bool MGSurfCurve::operator==(const MGSurfCurve& crv)const{
	return (m_surface==crv.m_surface) && (m_curve==crv.m_curve);
}
bool MGSurfCurve::operator<(const MGSurfCurve& gel2)const{
	if(m_surface==gel2.m_surface)
		return m_curve<gel2.m_curve;
    return (*m_surface)<*(gel2.m_surface);
}
bool MGSurfCurve::operator==(const MGGel& gel2)const{
	const MGSurfCurve* gel2_is_this=dynamic_cast<const MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGSurfCurve::operator<(const MGGel& gel2)const{
	const MGSurfCurve* gel2_is_this=dynamic_cast<const MGSurfCurve*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//Obtain parameter value if this curve is negated by "negate()".
double MGSurfCurve::negate_param(double t)const{
	assert(false);
	return 0.0;
}

MGPosition MGSurfCurve::negate_param(const MGPosition& t)const{
	assert(false);
	return MGPosition();
}

//Test if given point is on the curve or not. If yes, return parameter
//value of the curve. Even if not, return nearest point's parameter.
// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the curve,
// and 0 if the point is not on the curve.
/*bool MGSurfCurve::on(
	const MGPosition& P,//Point(指定点)
	double&	t			//Parameter of the curve(パラメータ)
	) const{
	double t0=param_s(), t1=param_e();
	MGVector V(eval(t0)-P);
	double lensqrMin0=V%V, tmin=t0;
	double lensqrMin=lensqrMin0;
	V=eval(t1)-P;
	double lensqr=V%V;
	if(lensqr<lensqrMin){
		lensqrMin=lensqr; tmin=t1;
	}

	size_t ndiv=intersect_dnum();
	double tg=t0;
	double delta=(t1-t0)/ndiv;

	for(size_t i=1; i<ndiv; i++){
		tg+=delta;
		V=eval(tg)-P;
		lensqr=V%V;
		if(lensqr<lensqrMin){lensqrMin=lensqr; tmin=tg;}
	}
	if(perp_guess(1.,0.,P,tmin,t)){
		V=eval(t)-P;
		lensqr=V%V;
		if(lensqr<=lensqrMin) lensqrMin=lensqr;
		else t=tmin;
	}else{
		//Try again by increasing ndiv.
		ndiv*=3;
		delta=(t1-t0)/ndiv;
		tg=t0;
		lensqrMin=lensqrMin0;
		for(i=1; i<ndiv; i++){
			tg+=delta;
			V=eval(tg)-P;
			lensqr=V%V;
			if(lensqr<lensqrMin){lensqrMin=lensqr; tmin=tg;}
		}
		if(perp_guess(1.,0.,P,tmin,t)){
			V=eval(t)-P;
			lensqr=V%V;
			if(lensqr<=lensqrMin) lensqrMin=lensqr;
			else t=tmin;
		}else t=tmin;
	}
	return lensqrMin<=MGTolerance::wc_zero_sqr();
}*/

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGSurfCurve::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGStraight& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGRLBRep& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGEllipse& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGLBRep& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGSurfCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}
MGPosition_list MGSurfCurve::perps(
	const MGBSumCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	perps_of_each(crv2,list);
	return list;	
}

//isect of this SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list MGSurfCurve::perps_noCompo(const MGCurve& curve2)const{
	MGPosition_list list=curve2.perps_with_noCompoSC(*this);
	return MGPosition_list(list,1,0);	
}

//perpendicular points of each elements of this m_curve,
//if m_curve is MGTrimmedCurve of a MGCompositeCurve.
void MGSurfCurve::perps_of_each(
	const MGCurve& curve2,	//The perps objective curve.
	MGPosition_list& list	//Obtained perpendicular points will be appended.
)const{
	const MGCompositeCurve* compo=base_composite();
	if(!compo){
		list=perps_noCompo(curve2);
		return;
	}

	MGCompositeCurve::const_iterator
		icurve=compo->begin(), endcurve=compo->end();
	double ts=m_curve.param_s(), te=m_curve.param_e();
	for(; icurve!=endcurve; endcurve++){
		const MGCurve& crvi=**icurve;
		if(crvi.param_e()<=ts || te<=crvi.param_s())
			continue;
		MGSurfCurve scrv(*m_surface,crvi);
		MGPosition_list listi=scrv.perps_noCompo(curve2);
		MGPosition_list::iterator j=listi.begin(), je=listi.end();
		for(; j!=je; j++){
			if(in_range(j->ref(0)))
				list.append(*j);
		}
	}
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGSurfCurve& MGSurfCurve::coordinate_exchange(size_t i, size_t j){
	assert(false);
	return *this;
}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGSurfCurve::copy_as_nurbs()const{return new MGLBRep(*this);}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGCurve* MGSurfCurve::copy_change_dimension(
	size_t sdim,			// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2		// Source order of this line.
)const{
	MGLBRep tmp(*this);
	return tmp.copy_change_dimension(sdim,start1,start2);
}

//Update this by limiting the parameter range of the curve.
// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
MGSurfCurve& MGSurfCurve::limit(const MGInterval& I){
	m_curve.limit(I);
	update_mark();
	return *this;
}

//Negate the curve direction(曲線の方向を反転する)
void MGSurfCurve::negate(){
	assert(false);
}

//Unlimit parameter range of the curve(limitをはずす)
MGCurve& MGSurfCurve::unlimit(){
	m_curve.unlimit();
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGSurfCurve::unlimit_end(){
	m_curve.unlimit_end();
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGSurfCurve::unlimit_start(){
	m_curve.unlimit_start();
	update_mark();
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線を平行移動して自身とする。
MGSurfCurve& MGSurfCurve::operator+= (const MGVector& v){
//	assert(false);
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGSurfCurve& MGSurfCurve::operator-= (const MGVector& v){
//	assert(false);
	return *this;
}

//Update the curve by multiplying scale.
// 与えられたスケールを曲線にかける。
MGSurfCurve& MGSurfCurve::operator*= (double scale){
//	assert(false);
	return *this;
}

//Update the curve by transformation of matrix.
// 与えられた変換で直線の変換を行い自身の直線とする。
MGSurfCurve& MGSurfCurve::operator*= (const MGMatrix& mat){
//	assert(false);
	return *this;
}

//Update the curve by transformation of transf.
// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGSurfCurve& MGSurfCurve::operator*= (const MGTransf& tr){
//	assert(false);
	return *this;
}

MGSurfCurve MGSurfCurve::operator+ (const MGVector& v) const{
	assert(false); return *this;
}
MGSurfCurve operator+ (const MGVector& v, const MGSurfCurve& cv2){
	assert(false); return cv2;
}
MGSurfCurve MGSurfCurve::operator- (const MGVector& v) const{
	assert(false); return *this;
}
MGSurfCurve MGSurfCurve::operator* (double scale) const{
	assert(false); return *this;
}
MGSurfCurve operator* (double scale, const MGSurfCurve& cv2){
	assert(false); return cv2;
}
MGSurfCurve MGSurfCurve::operator* (const MGMatrix& mat) const{
	assert(false); return *this;
}
MGSurfCurve MGSurfCurve::operator* (const MGTransf& tr) const{
	assert(false); return *this;
}
