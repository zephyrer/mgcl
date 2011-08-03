/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/Matrix.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Surface.h"
#include "mg/CParam_list.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/Plane.h"
#include "mg/SPhere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"
#include "topo/Face.h"

extern "C" {
#include "cskernel/Bleval.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGSBRep3.cpp
//
// Implements Surface B-Representation class MGSBRep.

//Intersection computation.

//The following two function will be used in perps or isect
//to decide how many division of the surface along u or v direction
//should be applied before using perp_guess or isect_guess.
size_t MGSBRep::intersect_dnum_u() const{
	return (bdim_u()+2-order_u())*(order_u()-1);
}
size_t MGSBRep::intersect_dnum_v() const{
	return (bdim_v()+2-order_v())*(order_v()-1);
}

// Surface と Curve の交点を求める。
MGCSisect_list MGSBRep::isect(const MGCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
MGCSisect_list MGSBRep::isect(const MGRLBRep& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
MGCSisect_list MGSBRep::isect(const MGEllipse& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
MGCSisect_list MGSBRep::isect(const MGLBRep& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
MGCSisect_list MGSBRep::isect(const MGSurfCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
MGCSisect_list MGSBRep::isect(const MGBSumCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Surface の交線を求める。
//Compute intesection of Sphere and Surface.
MGSSisect_list MGSBRep::isect(const MGSurface& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGSBRep::isect(const MGPlane& srf2) const{
	return intersectPl(srf2);
}
MGSSisect_list MGSBRep::isect(const MGSphere& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSBRep::isect(const MGCylinder& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSBRep::isect(const MGSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSBRep::isect(const MGRSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSBRep::isect(const MGBSumSurf& srf2) const{
	return intersect(srf2);
}

//Intersection of Surface and a straight line.
MGCSisect_list MGSBRep::isectSl(
	const MGStraight& sl,
	const MGBox& uvbox //indicates if this surface is restrictied to the parameter
					//range of uvbox. If uvbox.is_null(), no restriction.
)const{
	MGCSisect_list list(&sl,this);
	const MGBox& sbx=box();
	if(!sbx.crossing(sl))
		return list;

	MGUnit_vector SLD=sl.direction();
	MGMatrix mat; mat.to_axis(SLD,2);	//Matrix to transform SLD to be z axis.

	double u0,u1,v0,v1;
	bool uvbox_is_null=uvbox.is_null();
	if(uvbox_is_null){
		u0=param_s_u(); u1=param_e_u();
		v0=param_s_v(); v1=param_e_v();
	}else{
		u0=uvbox[0].low_point(); u1=uvbox[0].high_point();
		v0=uvbox[1].low_point(); v1=uvbox[1].high_point();
		double uspan=u1-u0, vspan=v1-v0;
		if(MGREqual_base(u0,param_s_u(),uspan) && 
			MGREqual_base(u1,param_e_u(),uspan) && 
			MGREqual_base(v0,param_s_v(),vspan) &&
			MGREqual_base(v1,param_e_v(),vspan) ) uvbox_is_null=true;
	}
	double u=(u0+u1)*.5, v=(v0+v1)*.5;
	MGVector PN=sl.eval(sl.closest(eval(u,v)));
		//PN is the nearest point of the sl from the center of this surface,
		//and is the vector to translate sl to pass through the origin.
	MGSBRep* sf2D=static_cast<MGSBRep*>(copy_surface());
	(*sf2D)-=PN; (*sf2D)*=mat;
	sf2D->change_dimension(2);
	//sf2D is the 2D surface that is transformed from the original this
	//surface so that the straight line sl is seen as a point of origin(0.,.0.).
	if(!uvbox_is_null) sf2D->limit(uvbox);
	//std::cout<<(*sf2D);////////////

	size_t m=sf2D->bdim_u(), n=sf2D->bdim_v();
	size_t mm1=m-1, nm1=n-1;
	MGSPointSeq& sp=sf2D->surface_bcoef();
	MGKnotVector& tu=sf2D->knot_vector_u();
	MGKnotVector& tv=sf2D->knot_vector_v();
	size_t kum1=tu.order()-1, kvm1=tv.order()-1;

	for(size_t i=kum1; i<m; i++){
	size_t ip1=i+1;
	while(tu(i)== tu(ip1) && ip1<m) {i=ip1; ip1++;}
	for(size_t j=kvm1; j<n; j++){
		size_t jp1=j+1;
		while(tv(j)== tv(jp1) && jp1<n) {j=jp1; jp1++;}
		MGBox bx;
		for(size_t iu=i-kum1; iu<=i; iu++)
			for(size_t ju=j-kvm1; ju<=j; ju++) bx.expand(sp(iu,ju));
		if(bx.includes_origin()){
			//Compute intersection point using isect_guess_straight.
			MGPosition uv((tu[i]+tu[ip1])*.5, (tv[j]+tv[jp1])*.5);
			double t=sl.closest(eval(uv));
			if(isect_guess_straight(sl,t,uv,t,uv)){
					if(uvbox_is_null || uvbox>>uv) list.append(sl.eval(t),t,uv);
			}
		}
	}
	}

	delete sf2D;
	return list;
}

//Return order of intersection line order of MGLBRep.
size_t MGSBRep::isect_order() const{
	size_t ku=order_u(), kv=order_v();
	if(ku<kv) ku=kv;
	return ku;
}

//"isect_sub_interval" is a dedicated function of isect_incr_pline.
//Computes subinterval id(index) and length of parameter line necessary for
//intersection computation(function's return value).
size_t MGSBRep::isect_sub_interval(
	size_t kdt,	//indicates if u=const(kdt=1) or v=const parameter line is used
				//or not to get the intersection
	double u, double v,	// u and v parameter value.
	double du, double dv,// incremental parameter length so far.
	size_t& index,		//index of B-Rep coef that is the start of
	//sub line b-rep(parameter line) for the computation of intersection
	size_t incr			//Incremental valuse of B-coef's id.
)const{
	double t0,dt;
	const MGKnotVector* tpointer;
	if(kdt){ //For u=const v-parameter line.
		t0=v; dt=dv; tpointer=&m_vknot;
	}else{	 //For v=const u-parameter line.
		t0=u; dt=du; tpointer=&m_uknot;
	}
	const MGKnotVector& t=*tpointer;

	int k=t.order(),n=t.bdim();
	int km1=k-1;
	int i, i1=t.locate(t0);
	if(incr>=1){
		if(dt>0.){
			i1+=incr; if(i1>=n) i1=n-1;
		}else{
			i1-=incr; if(i1<k) i1=km1;
		}
	}
	i=i1;
	while(i>=k && t(i)==t(i-1)) i--;
	int j=i1+3; 
	if(j>n) j=n; else while(j<n && t(j)==t(j-1)) j++;
	i=i-k-1; if(i<0) i=0;
	while(i>0 && t(i+km1)==t(i+k)) i--;
	int id=i;
	size_t cnum=j-i;
	index=id;
	return cnum;
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
void MGSBRep::isect_incr_pline2(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	size_t incr,			//Incremental valuse of B-coef's id.
	MGLBRep& pline			//parameter line will be output.
)const{
	const int nderiv=0, jcont=1;
	const size_t ncd=sdim(), ncd1=1;

	size_t i,j;
	size_t ku=order_u(), kv=order_v();
	size_t lud=bdim_u(), lvd=bdim_v();
	size_t isr1,isr2,isr12;
	m_surface_bcoef.capacity(isr1,isr2);isr12=isr1*isr2;
	double P_area[4]; double* P;
	if(ncd<=4) P=P_area; else P=new double[ncd];

	//Compute necessary sub-interval length of parameter line of this surface.
	size_t cnum,id;
	if(kdt){ u=m_uknot.range(uv[0]+du); v=uv[1];}
	else { v=m_vknot.range(uv[1]+dv); u=uv[0];}
	cnum=isect_sub_interval(kdt,u,v,du,dv,id,incr);

	const MGKnotVector* knot; if(kdt) knot=&m_vknot; else knot=&m_uknot;
	pline.knot_vector()=MGKnotVector(id,cnum,*knot);
	MGBPointSeq& lb=pline.line_bcoef();
	lb.resize(cnum,ncd);

	int ism=1, kar; 
	if(kdt){
	//Compute u=const v-parameter line of this surface in pline.
		kar=1;
		for(i=0; i<cnum; i++){
			bleval_(ku,lud,knot_data_u(),m_surface_bcoef.data(0,id+i,0),
				   isr12,ncd,kar,u,nderiv,jcont,ism,P);
			ism=2;
			for(j=0; j<ncd; j++) lb(i,j)=P[j];
		}
	}else{
	//Compute v=const u-parameter line in pline.
		kar=2;
		for(i=0; i<cnum; i++){
			for(j=0; j<ncd; j++){
				bleval_(kv,lvd,knot_data_v(),m_surface_bcoef.data(id+i,0,j),
					isr1,ncd1,kar,v,nderiv,jcont,ism,P+j);
				ism=2;
				lb(i,j)=P[j];
			}
		}
	}
	if(ncd>4) delete[] P;
}

//"isect_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGSBRep::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,	//Incremental parameter length.
	double& u, double& v,	//next u, v value will be output.
	size_t incr			//Invremental id value of B-coef's.
)const{
	MGLBRep* lbrep=new MGLBRep;
	isect_incr_pline2(uv,kdt,du,dv,u,v,incr,*lbrep);
	return lbrep;
}

//Obtain 1D surface rep. of this surf which can be used for
//isect(const MGPlane& pl). This surf1D is used in isect for
//the argument of isect_startPlane, which will use surf1D to compute isect(pl).
//surf1D=0.(intersection with x=0. plane) is the intersection lines.
MGSBRep* MGSBRep::surf1D(
	const MGPlane& pl
)const{
	MGVector cntr=pl.eval(pl.param(center()));
	MGSBRep* surf=new MGSBRep;
	//Construct a working 1D surface sf1D.
	MGUnit_vector N=MGVector(3,pl.normal());//Plane Normal(Unit-vector).
		//PN is the vector to translate pl to pass through the origin.
	MGSPointSeq cp=surface_bcoef()-cntr;
	MGMatrix mat; mat.to_axis(N,2);	//Matrix to transform N to be z axis.
	cp*=mat;

	surf->m_surface_bcoef=MGSPointSeq(1,cp,0,2);
	surf->m_uknot=knot_vector_u();
	surf->m_vknot=knot_vector_v();
	return surf;
}
