/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <stdlib.h>
#include "mg/Position_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Transf.h"
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

using namespace std;

//MGSphere is a Sphere in 3D space. Sphere f(u,v) is expressed
// by two ellipses EL1(m_ellipseu) and EL2(m_ellipsev) as:
//f(u,v) = C+(M*cos(u)+N*sin(u))*cos(v)+B*sin(v), or
//f(u,v) = C+EL1(u)*cos(v)+B*sin(v),
//where EL1=M*cos(u)+N*sin(u), and EL2=C+N*cos(v)+B*sin(v).
//Here M is the major axis of EL1, N is the minor axis of EL1, N is
//also the major axis of EL2, and B=(unit vector of (M*N))*(N.len()),
//which is the minor axis of EL2. (M,N,B) make a orthonormal system.
//v is the angle with M axis in the (B,M) plane, and u is the angle with M in the (M,N) plane.
//v=0 parameter line makes the ellipse C+EL1, and u=pai/2 parameter line
//makes the ellipse EL2.

// MGSphereクラスは３次元空間における球を表すクラスである。

/////////////Constructor コンストラクタ////////////

//Void constructor 初期化なしで柱面を生成する。
MGSphere::MGSphere(void):MGSurface(){;}

// Construct a  whole sphere from the center and the radius.
MGSphere::MGSphere(
	const MGPosition& cntr,	// Sphere center.
	double radius) 		// Sphere radius.
:MGSurface(){
	MGVector M=mgX_UVEC*radius;
	MGVector N=mgY_UVEC*radius;
	MGVector B=mgZ_UVEC*radius;
	m_ellipseu=MGEllipse(mgORIGIN,M,N,MGInterval(0.,mgDBLPAI));
	m_ellipsev=MGEllipse(cntr,N,B,MGInterval(-mgHALFPAI,mgHALFPAI));
}

// Construct a  whole sphere from the center and the radius.
//Let MGUnit_vector N(B*M), M2(N*B). Then (M2,N,B) makes a orthonormal system,
//and this sphere is parameterized as:
//F(u,v)=cntr+radis*cos(v)(M*cos(u)+N*sin(u))+radis*sin(v)*B.
MGSphere::MGSphere(
	const MGPosition& cntr,	// Sphere center.
	double radius,			// Sphere radius.
	const MGUnit_vector& B,	//axis
	const MGVector& M //reference direciotn that is approximately perpendiculat to B.
):MGSurface(){
	MGVector B2=B*radius;
	MGUnit_vector N2U=B*M; MGVector N2=N2U*radius;
	MGUnit_vector M2U=N2U*B; MGVector M2=M2U*radius;
	m_ellipseu=MGEllipse(mgORIGIN,M2,N2,MGInterval(0.,mgDBLPAI));
	m_ellipsev=MGEllipse(cntr,N2,B2,MGInterval(-mgHALFPAI,mgHALFPAI));
}

// Construct a Sphere by changing this space dimension or
// ordering the coordinates.
MGSphere::MGSphere(
	size_t dim,				// New space dimension.
	const MGSphere& sphr,	// Original Sphere.
	size_t start1, 		// Destination order of new Surface.
	size_t start2) 		// Source order of original Surface.
:MGSurface(sphr),
m_ellipseu(dim,sphr.m_ellipseu,start1,start2),
m_ellipsev(dim,sphr.m_ellipsev,start1,start2){
	update_mark();
}

//Whole sphere(u parameter range is from 0 to 2 pai) around minor axis of
//the input ellipse. The parameter range of the ellipse must be within the range
//from -pai/2 to pai/2, and the range becomes the v parameter range of the sphere.
//The input ellipse makes u=const(u=pai/2) v-paramarter line.
MGSphere::MGSphere(
	const MGEllipse& ellipse	//ellispe  of the Sphere
):MGSurface(),m_ellipsev(ellipse){
	m_ellipsev.limit(MGInterval(-mgHALFPAI,mgHALFPAI));
	const MGVector& N=m_ellipsev.major_axis();
	const MGVector& B=m_ellipsev.minor_axis();
	MGUnit_vector M=N*B;
	m_ellipseu=MGEllipse(mgORIGIN,M*(N.len()),N,MGInterval(0.,mgDBLPAI));
}

//Sphere(u parameter range is urange) around minor axis of the input ellipse.
//The parameter range of the ellipse must be in the range from -pai/2 to pai/2,
//and the range becomes the v parameter range of the sphere.
//The input ellipse makes u=const(u=pai/2) v-paramarter line.
MGSphere::MGSphere(
	const MGEllipse& ellipse,	//ellispe  of the Sphere
	MGInterval urange
):MGSurface(),m_ellipsev(ellipse){
	m_ellipsev.limit(MGInterval(-mgHALFPAI,mgHALFPAI));
	const MGVector& N=m_ellipsev.major_axis();
	const MGVector& B=m_ellipsev.minor_axis();
	MGUnit_vector M=N*B;
	m_ellipseu=MGEllipse(mgORIGIN,M*(N.len()),N,urange);
}

//////////Operator overload 演算子の多重定義/////////////

//Assignment.
//When the leaf object of this and srf2 are not equal, this assignment
//does nothing.
MGSphere& MGSphere::operator=(const MGSphere& gel2){
	if(this==&gel2)
		return *this;

	MGSurface::operator=(gel2);
	m_ellipseu=gel2.m_ellipseu;
	m_ellipsev=gel2.m_ellipsev;
	return *this;
}
MGSphere& MGSphere::operator=(const MGGel& gel2){
	const MGSphere* gel2_is_this=dynamic_cast<const MGSphere*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Translation of the Sphere
MGSphere MGSphere::operator+ (const MGVector& vec)const{
	MGSphere sphr(*this);
	sphr+=vec;
	return sphr;
}
MGSphere operator+ (const MGVector& v, const MGSphere& sphr){
	return sphr+v;
}

//Translation of the Sphere
MGSphere& MGSphere::operator+= (const MGVector& vec){
	m_ellipsev+=vec;
	if(m_box)
		(*m_box)+=vec;
	return *this;
}

//Translation of the Sphere
MGSphere MGSphere::operator- (const MGVector& vec) const{
	MGSphere sphr(*this);
	sphr-=vec;
	return sphr;
}

//Translation of the Sphere
MGSphere& MGSphere::operator-= (const MGVector& vec){
	m_ellipsev-=vec;
	if(m_box)
		(*m_box)-=vec;
	return *this;
}

//柱面のスケーリングを行い，柱面を作成する。
//Scaling of the Sphere by a double.
MGSphere MGSphere::operator* (double s) const{
	MGSphere sphr(*this);
	sphr*=s;
	return sphr;
}

//Scaling of the Sphere by a double.
MGSphere operator* (double scale, const MGSphere& sphr){
	return sphr*scale;
}

//Scaling of the Sphere by a double.
MGSphere& MGSphere::operator*= (double s){
	m_ellipsev*=s;
	m_ellipseu*=s;
	update_mark();
	return *this;
}

//Transformation of the Sphere by a matrix.
MGSphere MGSphere::operator* (const MGMatrix& mat) const{
	MGSphere sphr(*this);
	sphr*=mat;
	return sphr;
}

//Transformation of the Sphere by a matrix.
MGSphere& MGSphere::operator*= (const MGMatrix& mat){
	m_ellipseu*=mat;
	m_ellipsev*=mat;
	update_mark();
	return *this;
}

//Transformation of the Sphere by a MGTransf.
MGSphere MGSphere::operator* (const MGTransf& tr) const{
	MGSphere sphr(*this);
	sphr*=tr;
	return sphr;
}

//Transformation of the Sphere by a MGTransf.
MGSphere& MGSphere::operator*= (const MGTransf& tr){
	m_ellipseu*=tr.affine();
	m_ellipsev*=tr;
	update_mark();
	return *this;
}

//Comparison between Sphere and a surface.
bool MGSphere::operator==(const MGSphere& sphr)const{
	if(m_ellipsev!=sphr.m_ellipsev)
		return false;
	if(m_ellipseu!=sphr.m_ellipseu)
		return false;

	return true;
}
bool MGSphere::operator<(const MGSphere& gel2)const{
	if(m_ellipseu==gel2.m_ellipseu)
		return m_ellipsev<gel2.m_ellipsev;
	return m_ellipseu<gel2.m_ellipseu;
}
bool MGSphere::operator==(const MGGel& gel2)const{
	const MGSphere* gel2_is_this=dynamic_cast<const MGSphere*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGSphere::operator<(const MGGel& gel2)const{
	const MGSphere* gel2_is_this=dynamic_cast<const MGSphere*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

////////////Member function メンバ関数///////////

//Compute this surface's box
void MGSphere::box_driver(MGBox& bx)const{
	bx.set_null();

	box_vconst(bx,param_s_v());
	box_vconst(bx,param_e_v());
	box_uconst(bx,param_s_u());
	box_uconst(bx,param_e_u());

	const MGVector& m=M();
	const MGVector& n=N();

	for(size_t i=0; i<3; i++){

	double ni=n[i], mi=m[i];
	double sqrtmini=sqrt(mi*mi+ni*ni);
	if(MGMZero(sqrtmini))
		continue;

	double ani=fabs(ni), ami=fabs(mi);
	double ui;
	if(ani<=ami)
		ui=asin(ani/sqrtmini);
	else
		ui=acos(ami/sqrtmini);

	if(m_ellipseu.in_RelativeRange_of_radian(ui)){
		double u=m_ellipseu.RelativeRange_in_radian(ui);
		box_uconst(bx,u);
	}
	if(m_ellipseu.in_RelativeRange_of_radian(-ui)){
		double u=m_ellipseu.RelativeRange_in_radian(-ui);
		box_uconst(bx,u);
	}
	if(m_ellipseu.in_RelativeRange_of_radian(ui-mgPAI)){
		double u=m_ellipseu.RelativeRange_in_radian(ui-mgPAI);
		box_uconst(bx,u);
	}
	if(m_ellipseu.in_RelativeRange_of_radian(mgPAI-ui)){
		double u=m_ellipseu.RelativeRange_in_radian(mgPAI-ui);
		box_uconst(bx,u);
	}

	}
}

// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
//Box that includes limitted Sphere by box.
MGBox MGSphere::box_limitted(
	const MGBox& uvrange	// Parameter Range of the surface.
) const{
	MGSphere sphr(*this);
	sphr.m_ellipseu.limit(uvrange[0]);
	sphr.m_ellipsev.limit(uvrange[1]);
	return sphr.box();
}

//Compute box of u=const parameter line.
//The box will be or-ed to the input box bx.
//u must be the value in radian.
void MGSphere::box_uconst(MGBox& bx, double u)const{
	MGVector MN(m_ellipseu.eval(u));
	double v0=m_ellipsev.gp_to_radian(m_ellipsev.param_s());
	double v1=m_ellipsev.gp_to_radian(m_ellipsev.param_e());
	MGEllipse elv(C(),MN,B(),MGInterval(v0,v1));
	bx|=elv.box();
}

//Compute box of v=const parameter line.
//The box will be or-ed to the input box bx.
//v must be the value in radian.
void MGSphere::box_vconst(MGBox& bx, double v)const{
	double vR=m_ellipsev.gp_to_radian(v);
	double cosv=cos(vR), sinv=sin(vR);
	MGEllipse elu(m_ellipseu*cosv+(C()+B()*sinv));
	bx|=elu.box();
}

//Changing this object's space dimension.
MGSphere& MGSphere::change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	m_ellipseu.change_dimension(sdim,start1,start2);
	m_ellipsev.change_dimension(sdim,start1,start2);
	update_mark();
	return *this;
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
MGSphere& MGSphere::change_range(
	int is_u,				//if true, (t1,t2) are u-value. if not, v.
	double t1,				//Parameter value for the start of original. 
	double t2				//Parameter value for the end of original. 
){
	if(is_u)
		m_ellipseu.change_range(t1,t2);
	else
		m_ellipsev.change_range(t1,t2);
	return *this;
}

//Compute the closest point parameter value (u,v) of this surface
//from a point.
MGPosition MGSphere::closest(const MGPosition& point) const{
	MGPosition_list list=perps(point);
	size_t n=list.size();
	if(n==2) return list.front();
	else if(n==1){
		MGPosition& Q=list.front();
		const MGPosition& C=sphere_center();
		if((Q-C)%(point-C)>0.) return Q;
	}

	//Compute using closest_on_perimeter().
	return closest_on_perimeter(point);
}

//Compute the closest point on a perimeter of the surface. The point is returned
//as the parameter value (u,v) of this surface.
MGPosition MGSphere::closest_on_perimeter(const MGPosition& point)const{
	MGPosition uv1, uv;
	double dist1,dist;
	size_t pnum=perimeter_num();
	MGCurve* perimtr;
	for(size_t i=0; i<pnum; i++){
		if(i==0 && degenerate_at_v0()) continue;
		if(i==2 && degenerate_at_v1()) continue;
		perimtr=perimeter_curve(i);
		uv1=perimeter_uv(i,perimtr->closest(point));
		dist1=(point-eval(uv1)).len();
		if(uv.is_null()){dist=dist1; uv=uv1;}
		else if(dist1<dist){dist=dist1; uv=uv1;}
		delete perimtr;
	}
	return uv;
}

//Return minimum box that includes whole of the surface.
//Returned is a newed object pointer.
MGBox* MGSphere::compute_box() const{
	MGBox* bx=new MGBox();
	box_driver(*bx);
	return bx;
}

//Compute the sphere parameter vaule (u,v) whose world coordinate is P,
//given a point P on the surface.
//Computed (u,v) may be outside real parameter range
//(sphere is regarded as a whole one).
void MGSphere::compute_uv(const MGPosition& P, double&u, double&v)const{
	assert(sphere());//******Currently this is valid only for sphere*******
	MGVector CP(P,sphere_center());//Vector from the center to P.
	const MGVector& Ms=M();
	const MGVector& Ns=N();
	const MGVector& Bs=B();
	v=mgHALFPAI-CP.angle(Bs);
	if(CP.orthogonal(Ms)){
		u=mgHALFPAI;
	}else{
		MGVector CPM=Bs*CP*Bs;//CPM is the projection vector of CP onto (M, N) plane.
		u=CPM.angle2pai(Ms,Bs);
	}
	u=m_ellipseu.radian_to_gp(u);
	v=m_ellipsev.radian_to_gp(v);
}

//Construct new surface object by copying to newed area.
//User must delete this copied object by "delete".
MGSphere* MGSphere::clone() const{
	return new MGSphere(*this);
}

//Construct new surface object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGSphere* MGSphere::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2) 		// Source order of this line.
const{
	return new MGSphere(sdim,*this,start1,start2);
}

//Ask if this sphere has the degenerate point at v=min, or v=max.
bool MGSphere::degenerate_at_v0()const{//at v=min
	double v0=m_ellipsev.gp_to_radian(m_ellipsev.param_s());
	return MGREqual_base(v0,-mgHALFPAI,mgPAI);
}
bool MGSphere::degenerate_at_v1()const{//at v=max
	double v1=m_ellipsev.gp_to_radian(m_ellipsev.param_e());
	return MGREqual_base(v1,mgHALFPAI,mgPAI);
}
// 与えられた点との距離を返却する。
//Return the distace between Sphere and the point.
double MGSphere::distance(const MGPosition& point) const{
	return (eval(closest(point))-point).len();
}

//Evaluate surface data.
MGVector MGSphere::eval(
	double u, double v		// Parameter value of the surface.
	, size_t ndu			// Order of derivative along u.
	, size_t ndv			// Order of derivative along v.
)const{
	MGVector Fu=m_ellipseu.eval(u,ndu);
	double vr=m_ellipsev.gp_to_radian(v);
	double cosv=cos(vr), sinv=sin(vr);
	const MGVector& Bs=B();
	if(ndu==0 && ndv==0)
		return C()+Fu*cosv+Bs*sinv;

	MGVector F;
	div_t result=div(ndv, 4);
	if(ndu==0){
		switch(result.rem){
		case 0: F=  Fu*cosv+Bs*sinv; break;
		case 1: F= -Fu*sinv+Bs*cosv; break;
		case 2: F= -Fu*cosv-Bs*sinv; break;
		case 3: F=  Fu*sinv-Bs*cosv; break;
		}
	}else{
		switch(result.rem){
		case 0: F=  Fu*cosv; break;
		case 1: F= -Fu*sinv; break;
		case 2: F= -Fu*cosv; break;
		case 3: F=  Fu*sinv; break;
		}
	}
	return F;
}

// Exchange parameter u and v.
//This is not allowed.
MGSurface& MGSphere::exchange_uv(){ assert(false);return *this;}

//Modify the original Surface by extrapolating the specified perimeter.
//The extrapolation is C2 continuous if the order >=4.
//The extrapolation is done so that extrapolating length is "length"
//at the position of the parameter value "param" of the perimeter.
MGSphere& MGSphere::extend(
	int perimeter,	//perimeter number of the Surface.
					// =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	// parameter value of above perimeter.
	double length,	//chord length to extend at the parameter param of the perimeter.
	double dk){  //Coefficient of how curvature should vary at
//    extrapolation start point. When dk=0, curvature keeps same, i.e.
//    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
//    i.e. dK/dS=-K/length at extrapolation start point.
//    (S=parameter of arc length, K=Curvature at start point)
//    That is, when dk reaches to 1 from 0, curve changes to flat.

	assert(0<=perimeter && perimeter<4);
	
	bool ise=true;//starting perimeter
	MGEllipse* el_to_extend;
	size_t ndu=0,ndv=0;
	if(perimeter==1 || perimeter==3){	// Extrapolate to u-direction
		el_to_extend=&m_ellipseu;
		if(perimeter==1)
			ise=false;//ending perimeter
		ndu=1;
	}else{
		// Extrapolate to v-direction
		el_to_extend=&m_ellipsev;
		if(perimeter==2)
			ise=false;//ending perimeter
		ndv=1;
	}
	MGPosition uv=perimeter_uv(perimeter,param);//Surface parameter value of param.
	double vlen=eval(uv,ndu,ndv).len();
	if(MGMZero(vlen))
		return *this;

	double slen=length/vlen;
	el_to_extend->extend(slen,ise);
	if(ndv)
		el_to_extend->limit(MGInterval(-mgHALFPAI, mgHALFPAI));
		//When extension is done about m_ellisev, the parameter range is limitted.

	return *this;
}

bool MGSphere::in_range(double u, double v) const{
	return m_ellipseu.in_range(u) && m_ellipsev.in_range(v);
}

// Surface と Curve の交点を求める。
//Compute intesection of Sphere and Curve.
MGCSisect_list MGSphere::isect(const MGCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Sphere and Curve.
MGCSisect_list MGSphere::isect(const MGRLBRep& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Sphere and Curve.
MGCSisect_list MGSphere::isect(const MGEllipse& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Sphere and Curve.
MGCSisect_list MGSphere::isect(const MGLBRep& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Sphere and Curve.
MGCSisect_list MGSphere::isect(const MGSurfCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Sphere and Curve.
MGCSisect_list MGSphere::isect(const MGBSumCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Surface の交線を求める。
//Compute intesection of Sphere and Surface.
MGSSisect_list MGSphere::isect(const MGSurface& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGSphere::isect(const MGPlane& srf2) const{
	return intersectPl(srf2);
}
MGSSisect_list MGSphere::isect(const MGSphere& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSphere::isect(const MGCylinder& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSphere::isect(const MGSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSphere::isect(const MGRSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGSphere::isect(const MGBSumSurf& srf2) const{
	return intersect(srf2);
}

#define MGSphere_isect_div_num 12.
//isect_dt computes incremental values of u and v direction for the intersection
//computation at parameter position (u,v).
void MGSphere::isect_dt(
	double u, double v, double& du, double& dv,
	double acuRatio		//acuracy ratio.
)const{
	double u0=m_ellipseu.param_s(); u0=m_ellipseu.gp_to_radian(u0);
	double u1=m_ellipseu.param_e(); u1=m_ellipseu.gp_to_radian(u1);
	double v0=m_ellipsev.param_s(); v0=m_ellipsev.gp_to_radian(v0);
	double v1=m_ellipsev.param_e(); v1=m_ellipsev.gp_to_radian(v1);
	double alfa=isect_dt_coef(0);
	double ualfa=alfa*(mgDBLPAI/(u1-u0))/MGSphere_isect_div_num;
	double valfa=alfa*(mgDBLPAI/(v1-v0))/MGSphere_isect_div_num;
	du=m_ellipseu.param_span()*ualfa*acuRatio;
	dv=m_ellipsev.param_span()*valfa*acuRatio;
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGSphere::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	size_t incr		//Incremental valuse of B-coef's id.
) const{
	//Compute necessary sub-interval length of parameter line of this surface.
	if(kdt){
		u=m_ellipseu.range(uv[0]+du); v=uv[1];
		//Compute u=const v-parameter line of this surface in pline.
		return parameter_curve(kdt,u);
	}else{
		v=m_ellipsev.range(uv[1]+dv); u=uv[0];
		//Compute v=const u-parameter line in pline.
		return parameter_curve(kdt,v);
	}
}

//Compute the intersection line of this and the plane pl.
MGSSisect_list MGSphere::intersectPl(const MGPlane& pl)const{
	assert(sphere());//******Currently this is valid only for sphere*******
	MGSSisect_list list(this,&pl);
	if(fabs(pl.distance(C()))>radius()) return list;

	MGPosition_list uvuv_list;
	intersect12Boundary(pl,uvuv_list);
	if(!uvuv_list.size()) uvuv_list=intersectInner(pl);
	//Compute intersection line using isect_with_surf.
	return isect_with_surf(uvuv_list,pl);
}

//Intersection of Surface and a straight line.
MGCSisect_list MGSphere::isectSl(
	const MGStraight& sl,
	const MGBox& uvbox //indicates if this surface is restrictied to the parameter
					//range of uvbox. If uvbox.is_null(), no restriction.
)const{
	assert(sphere());//******Currently this is valid only for sphere*******
	MGCSisect_list list(&sl,this);

	double r=radius();
	if(sl.distance(C())>r) return list;
	if(sl.straight_type()==MGSTRAIGHT_SEGMENT){
		if((sl.start_point()-C()).len()<r && (sl.end_point()-C()).len()<r)return list;
	}
	MGStraight sl2(sl);sl2.unlimit();
	double tsl;
	sl2.perp_point(C(),tsl);
	//tsl=parameter of sl where the center of the sphere is perpendicular to the sl.
	MGPosition Q=sl2.eval(tsl);
	double d=(Q-C()).len();
	double l=sqrt(r*r-d*d);
	MGUnit_vector sldir(sl2.direction());
	MGPosition P1=Q+sldir*l, P2=Q-sldir*l;
		//Two intersection points with the infinite straight line .
	double t1,t2, u,v;
	if(sl.on(P1,t1)){
		compute_uv(P1,u,v);
		if(in_range(u,v)) list.append(P1,t1,MGPosition(u,v));
	}
	if(sl.on(P2,t2)){
		compute_uv(P2,u,v);
		if(in_range(u,v)) list.append(P2,t2,MGPosition(u,v));
	}
	return list;
}

//Negate the normal of the Sphere.
void MGSphere::negate(
	int is_u)// Negate along u-direction if is_u is ture,
			// else along v-direction.
{
	if(is_u)
		m_ellipseu.negate();
	else
		m_ellipsev.negate();
}

//Obtain parameter value if this surface is negated by "negate()".
//Negate along u-direction if is_u is ture,
// else along v-direction.
MGPosition MGSphere::negate_param(const MGPosition& uv, int is_u)const{
	MGPosition uvnew(uv);
	if(is_u)
		uvnew(0)=m_ellipseu.negate_param(uv[0]);
	else
		uvnew(1)=m_ellipsev.negate_param(uv[1]);
	return uvnew;
}

//C1連続曲面の一定オフセット関数
//オフセット方向は、ノーマル方向を正とする。トレランスはline_zero()を使用している。
//戻り値は、オフセット面のオートポインターが返却される。
//Surface offset. positive offset value is offset normal direction.
//the radius of curvature must be larger than offset value.
//line_zero() is used.
std::auto_ptr<MGSurface> MGSphere::offset_c1(
	double ofs_value,	//オフセット量
	int& error//エラーコード 0:成功 -1:面におれがある*
			// -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
)const{
	assert(sphere());//******Currently this is valid only for sphere*******

	if(!outgoing())
		ofs_value*=-1.;

	double r=radius();
	double scl=(r+ofs_value)/r;
	std::auto_ptr<MGSphere> surf(new MGSphere(*this));
	surf->m_ellipseu*=scl;
	surf->m_ellipsev*=scl;
	error=0;
	return surf;
}

// 指定点が面上にあるか調べる。（面上ならばtrue）
//Test if a point is on the Sphere. If on the Sphere, return true.
bool MGSphere::on(
	const MGPosition& point,	//A point to test 指定点
	MGPosition& puv				//Parameter value of the Sphere will be
								//returned.
)const{
	puv=closest(point);
	return ((point-eval(puv)).len()<=MGTolerance::wc_zero());
}

// Output virtual function.
//Output to ostream メンバデータを標準出力に出力する。
std::ostream& MGSphere::out(std::ostream &outpt) const{
	outpt<<"MGSphere::"<<this<<",";
	MGSurface::out(outpt);
	outpt<<std::endl<<" ,m_ellipseu="<<m_ellipseu<<std::endl;
	outpt<<" ,m_ellipsev="<<m_ellipsev;
	return outpt;
}

//test if the surface normal is outgoing from the center or not.
//If the sphere normao is outgoing, retrun true.
bool MGSphere::outgoing()const{
	MGPosition uvmid=param_mid();
	MGVector nrml=normal(uvmid);
	MGVector CP=eval(uvmid)-C();// the vector from center to mid point.
	return CP%nrml>0.;
}

//Obtain parameter space error.
double MGSphere::param_error() const{
	double uerror=m_ellipseu.param_error();
	double verror=m_ellipsev.param_error();
	return sqrt(uerror*uerror+verror*verror);
}

// パラメータ範囲を返す。
//Return parameter range of the Sphere(Infinite box).
MGBox MGSphere::param_range() const{
	return MGBox(m_ellipseu.param_range(), m_ellipsev.param_range());
}

// Compute parameter curve.
//Returned is newed area pointer, and must be freed by delete.
MGCurve* MGSphere::parameter_curve(
	int is_u	//Indicates x is u-value if is_u is true.
	, double x	//Parameter value.
				//The value is u or v according to is_u.
)const{
	if(is_u){
		MGVector MN(m_ellipseu.eval(x));
		double v0=m_ellipsev.param_s();
		double v1=m_ellipsev.param_e();
		double v0R=m_ellipsev.gp_to_radian(v0);
		double v1R=m_ellipsev.gp_to_radian(v1);
		MGEllipse* el=new MGEllipse(C(),MN,B(),MGInterval(v0R,v1R));
		el->change_range(v0,v1);
		return el;
	}else{
		double vR=m_ellipsev.gp_to_radian(x);
		double cosv=cos(vR), sinv=sin(vR);
		MGEllipse* el=new MGEllipse(m_ellipseu*cosv+(C()+B()*sinv));
		return el;
	}
}

//Compute part of the surface limitted by the parameter range bx.
//bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
MGSphere* MGSphere::part(
	const MGBox& uvbx,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
)const{
	MGSphere* sphr=new MGSphere(*this);
	sphr->m_ellipseu.limit(uvbx[0]);
	sphr->m_ellipsev.limit(uvbx[1]);
	return sphr;
}

// i must be < perimeter_num().
//When perimeter_num()==0, this function is undefined.
MGCurve* MGSphere::perimeter_curve(size_t i) const{
	if(i==0){
		return parameter_curve(0,param_s_v());
	}else if(i==2){
		return parameter_curve(0,param_e_v());
	}else if(i==1){
		return parameter_curve(1,param_e_u());
	}else{
		return parameter_curve(1,param_s_u());
	}
}

// 与えられた点にもっとも近い面上の垂直のパラメータ値を返却する。
//Return the nearest perpendicular point of the Sphere from P.
// Function's return value is whether point is obtained(1) or not(0)
int MGSphere::perp_point(
	const MGPosition& P,// 与えられた点
	MGPosition& uv,		//Parameter value of the Sphere will be output
	const MGPosition* uvguess	// guess
)const{
	MGPosition_list list=perps(P);
	if(list.size()){
		uv=list.front();//1st one is the nearer point from P.
		return 1;
	}
	return 0;
}

//Return all(actually at most two) foots of perpendicular straight lines from P.
//When two points are output, the nearer point will be output first
//in MGPosition_list.
MGPosition_list MGSphere::perps(
	const MGPosition& point				// Point in a space(指定点)
) const{
	assert(sphere());
	MGPosition_list list;
	const MGPosition& cntr=sphere_center();
	MGVector dir=point-cntr;
	double len=dir.len();
	if(MGMZero(len)){//when point is on the center of the sphere.
		list.append(MGPosition(param_s_u(), param_s_v()));
		list.append(MGPosition(param_e_u(), param_e_v()));
		return list;
	}
	double r=radius();
	dir*=r/len;
	MGPosition P1=cntr+dir, P2=cntr-dir;//The point on the sphere.
	double u,v;
	compute_uv(P1,u,v); if(in_range(u,v)) list.append(MGPosition(u,v));
	compute_uv(P2,u,v); if(in_range(u,v)) list.append(MGPosition(u,v));
	return list;
}

// 入力パラメータをパラメータ範囲でまるめて返却する。
//Round the input uv into parameter range of the Sphere, 
//return the same value as input.
MGPosition MGSphere::range(const MGPosition& uv) const{
	return MGPosition(m_ellipseu.range(uv[0]), m_ellipsev.range(uv[1]));
}

//Return the space dimension.
size_t MGSphere::sdim() const{return 3;}

//メンバデータを読み込む関数
void MGSphere::ReadMembers(MGIfstream& buf){
	MGSurface::ReadMembers(buf);
	m_ellipsev.ReadMembers(buf);
	m_ellipseu.ReadMembers(buf);
}

//メンバデータを書き込む関数
void MGSphere::WriteMembers(MGOfstream& buf) const{
	MGSurface::WriteMembers(buf);
	m_ellipsev.WriteMembers(buf);
	m_ellipseu.WriteMembers(buf);
}
