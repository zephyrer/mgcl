/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Transf.h"
#include "mg/Plane.h"
#include "mg/SPhere.h"
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

// MGCylinder.cpp
// Implementation of class MGCylinder
//MGCylinder is a Cylinder in 3D space.
//Cylinder is expressed by an ellipse and a straight line.
//Cylinder function  f(u,v) = m_ellipse(u) + m_axis(v),
//where u and v are two parameter of surface representation.
//Here, m_axis is a straight line that passes through the origin.
//m_axis'es m_root_point is always set to the origin.
//m_ellise is the ellipse when v=0(when v=0, m_axis(v)=the origin);

// MGCylinderクラスは３次元空間における円筒面を表すクラスである。
// MGCylinderクラスでは以下のようなパラメータ表現を使用します。
// f(u,v) = m_ellipse(u) + m_axis(v);

/////////////Constructor コンストラクタ////////////

//Void constructor 初期化なしで柱面を生成する。
MGCylinder::MGCylinder(void):MGSurface(){;}

// Construct a Cylinder by changing this space dimension or
// ordering the coordinates.
MGCylinder::MGCylinder(
	size_t dim,				// New space dimension.
	const MGCylinder& cyl,	// Original Cylinder.
	size_t start1, 		// Destination order of new Surface.
	size_t start2) 		// Source order of original Surface.
:MGSurface(cyl),m_ellipse(dim,cyl.m_ellipse,start1,start2),
m_axis(dim,cyl.m_axis,start1,start2){
	update_mark();
	m_ortho=m_ellipse.normal().parallel(m_axis.direction());
}

// Construct a cylinder of whole circle whose bottom center is bottom
//and top center is bottom+axis.
MGCylinder::MGCylinder(
	const MGPosition& bottom,	//Location on axis to position the cylinder,
						//defines zero v. 
	const MGVector axis,//The axis vector for the cylinder. 
	double radius,		//The radius of the cylinder.
	bool outgoing		//Indicates if the surface normal is going to
						//outside of the cylinder(true) or inside(false).
):MGSurface(),m_ellipse(bottom,radius,axis),
m_axis(mgORIGIN+axis,mgORIGIN), m_ortho(true){
	if(!outgoing) m_ellipse.negate();
}

//Cylinder from an ellipse and the axis straight line.
//When axis'es start point is not ellipse's center, axis'es start point
//will be set to the center.
MGCylinder::MGCylinder(
	const MGEllipse& ellipse,	//ellispe  of the cylinder
	const MGStraight& axis		//axis of the cylinder.
			//axis's root point will be neglected, always be set as
			//the origin.
):MGSurface(),m_ellipse(ellipse),m_axis(axis){
	m_axis.m_root_point=mgORIGIN;
	m_ortho=m_ellipse.normal().parallel(axis.direction());
}

//////////Operator overload 演算子の多重定義/////////////

//Assignment.
//When the leaf object of this and srf2 are not equal, this assignment
//does nothing.
MGCylinder& MGCylinder::operator=(const MGCylinder& gel2){
	if(this==&gel2)
		return *this;

	MGSurface::operator=(gel2);
	m_axis=gel2.m_axis;
	m_ellipse=gel2.m_ellipse;
	m_ortho=gel2.m_ortho;
	return *this;
}
MGCylinder& MGCylinder::operator=(const MGGel& gel2){
	const MGCylinder* gel2_is_this=dynamic_cast<const MGCylinder*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Translation of the Cylinder
MGCylinder MGCylinder::operator+ (const MGVector& vec) const{
	return MGCylinder(m_ellipse+vec,m_axis);
}
MGCylinder operator+ (const MGVector& v, const MGCylinder& cyl){
	return cyl+v;
}

//Translation of the Cylinder
MGCylinder& MGCylinder::operator+= (const MGVector& vec){
	m_ellipse+=vec;
	if(m_box) (*m_box)+=vec;
	return *this;
}

//Translation of the Cylinder
MGCylinder MGCylinder::operator- (const MGVector& vec) const{
	return MGCylinder(m_ellipse-vec,m_axis);
}

//Translation of the Cylinder
MGCylinder& MGCylinder::operator-= (const MGVector& vec){
	m_ellipse-=vec;
	if(m_box) (*m_box)-=vec;
	return *this;
}

//柱面のスケーリングを行い，柱面を作成する。
//Scaling of the Cylinder by a double.
MGCylinder MGCylinder::operator* (double s) const{
	return MGCylinder(m_ellipse*s, m_axis*s);
}

//Scaling of the Cylinder by a double.
MGCylinder operator* (double scale, const MGCylinder& cyl){
	return cyl*scale;
}

//Scaling of the Cylinder by a double.
MGCylinder& MGCylinder::operator*= (double s){
	m_ellipse*=s;
	m_axis*=s;
	update_mark();
	return *this;
}

//Transformation of the Cylinder by a matrix.
MGCylinder MGCylinder::operator* (const MGMatrix& mat) const{
	return MGCylinder(m_ellipse*mat, m_axis*mat);
}

//Transformation of the Cylinder by a matrix.
MGCylinder& MGCylinder::operator*= (const MGMatrix& mat){
	m_ellipse*=mat;
	m_axis*=mat;
	m_ortho=m_ellipse.normal().parallel(m_axis.direction());
	update_mark();
	return *this;
}

//Transformation of the Cylinder by a MGTransf.
MGCylinder MGCylinder::operator* (const MGTransf& tr) const{
	return MGCylinder(m_ellipse*tr, m_axis*tr.affine());
}

//Transformation of the Cylinder by a MGTransf.
MGCylinder& MGCylinder::operator*= (const MGTransf& tr){
	m_ellipse*=tr;
	m_axis*=tr.affine();
	m_ortho=m_ellipse.normal().parallel(m_axis.direction());
	update_mark();
	return *this;
}

//Comparison between Cylinder and a surface.
bool MGCylinder::operator==(const MGCylinder& srf2)const{
	if(m_ellipse!=srf2.m_ellipse)
		return false;
	if(m_axis!=srf2.m_axis)
		return false;

	return true;
}
bool MGCylinder::operator<(const MGCylinder& gel2)const{
	if(m_ellipse==gel2.m_ellipse)
		return m_axis<gel2.m_axis;
	return m_ellipse<gel2.m_ellipse;
}
bool MGCylinder::operator==(const MGGel& gel2)const{
	const MGCylinder* gel2_is_this=dynamic_cast<const MGCylinder*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGCylinder::operator<(const MGGel& gel2)const{
	const MGCylinder* gel2_is_this=dynamic_cast<const MGCylinder*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

/////////Debug function デバッグ関数///////////
// Output virtual function.
//Output to ostream メンバデータを標準出力に出力する。
std::ostream& MGCylinder::out(std::ostream &outpt) const{
	outpt<<"MGCylinder::"<<this<<",m_ortho="<<m_ortho<<std::endl<<" ,m_ellipes="<<m_ellipse
		<<" ,m_axis="<<m_axis;
	return outpt;
}

////////////Member function メンバ関数///////////

//Compute the axis point of the parameter v.
MGPosition MGCylinder::axis_point(double v)const{
	return m_ellipse.center()+m_axis.eval(v);
}

// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
//Box that includes limitted Cylinder by box.
MGBox MGCylinder::box_limitted(
	const MGBox& uvrange	// Parameter Range of the surface.
) const{
	MGEllipse el(m_ellipse); el.limit(uvrange[0]);
	MGStraight axis(m_axis); axis.limit(uvrange[1]);
	MGBox bx=el.box()+axis.start_point();
	bx|=(el+=axis.end_point()).box();
	return bx;
}

//Obtain ceter coordinate of the geometry.
MGPosition MGCylinder::center() const{
	return m_ellipse.center()+m_axis.center();
}

//Changing this object's space dimension.
MGCylinder& MGCylinder::change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	m_ellipse.change_dimension(sdim,start1,start2);
	m_axis.change_dimension(sdim,start1,start2);
	update_mark();
	return *this;
}

//Compute the closest point parameter value (u,v) of this surface
//from a point.
MGPosition MGCylinder::closest(const MGPosition& point) const{
	MGPosition uv;
	perp_point(point,uv);
	return uv;
}

//Return minimum box that includes whole of the surface.
//Returned is a newed object pointer.
MGBox* MGCylinder::compute_box() const{
	MGBox* bx=new MGBox(m_ellipse.box()+m_axis.start_point());
	*bx|=(m_ellipse+m_axis.end_point()).box();
	return bx;
}

//Construct new surface object by copying to newed area.
//User must delete this copied object by "delete".
MGCylinder* MGCylinder::clone() const{
	return new MGCylinder(*this);
}

//Construct new surface object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGCylinder* MGCylinder::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2) 		// Source order of this line.
const{
	return new MGCylinder(sdim,*this,start1,start2);
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
MGCylinder& MGCylinder::change_range(
	int is_u,				//if true, (t1,t2) are u-value. if not, v.
	double t1,				//Parameter value for the start of original. 
	double t2				//Parameter value for the end of original. 
){
	if(is_u)
		m_ellipse.change_range(t1,t2);
	else
		m_axis.change_range(t1,t2);
	return *this;
}

// 与えられた点との距離を返却する。
//Return the distace between Cylinder and the point.
double MGCylinder::distance(const MGPosition& point) const{
	return (eval(closest(point))-point).len();
}

//Evaluate surface data.
MGVector MGCylinder::eval(
	double u, double v		// Parameter value of the surface.
	, size_t ndu			// Order of derivative along u.
	, size_t ndv			// Order of derivative along v.
)const{
	if(ndu==0){
		if(ndv==0) return m_ellipse.eval(u)+m_axis.eval(v);
		else return m_axis.eval(v,ndv);
	}else{
		if(ndv==0) return m_ellipse.eval(u,ndu);
		else return mgZERO_VEC;
	}
}

//Evaluate right continuous surface data.
//Evaluate all positional data, 1st and 2nd derivatives.
void MGCylinder::eval_all(
	double u, double v,		// Parameter value of the surface.
	MGPosition& f,			// Positional data.
	MGVector&   fu,			// df(u,v)/du
	MGVector&   fv,			// df/dv
	MGVector&   fuv,		// d**2f/(du*dv)
	MGVector&   fuu,		// d**2f/(du**2)
	MGVector&   fvv			// d**2f/(dv**2)
)const{
	f=m_ellipse.eval(u)+m_axis.eval(v);
	fu=m_ellipse.eval(u,1);
	fv=m_axis.eval(v,1);
	fuv=fuu=fvv=mgZERO_VEC;
}

// Exchange parameter u and v.
//This is not allowed.
MGSurface& MGCylinder::exchange_uv(){ assert(false);return *this;}

//Modify the original Surface by extrapolating the specified perimeter.
//The extrapolation is C2 continuous if the order >=4.
//The extrapolation is done so that extrapolating length is "length"
//at the position of the parameter value "param" of the perimeter.
MGCylinder& MGCylinder::extend(
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
	
	bool is_start=true;//starting perimeter
	MGCurve* to_extend;
	size_t ndu=0,ndv=0;
	if(perimeter==1 || perimeter==3){	// Extrapolate to u-direction
		to_extend=&m_ellipse;
		if(perimeter==1)
			is_start=false;//ending perimeter
		ndu=1;
	}else{
		// Extrapolate to v-direction
		to_extend=&m_axis;
		if(perimeter==2)
			is_start=false;//ending perimeter
		ndv=1;
	}
	MGPosition uv=perimeter_uv(perimeter,param);//Surface parameter value of param.
	double vlen=eval(uv,ndu,ndv).len();
	double slen=length/vlen;
	to_extend->extend(slen,is_start);

	return *this;
}

//Test if the surface is flat or not within the parameter value rectangle of uvbox.
//Function's return value is:
//	true: if the surface is flat,
//  false: if the surface is not falt.
//When this is not falt, the direction that indicates which direction the surface
//should be divided will be output.
//***** the flatness is tested only approximately. This is for exclusive use of
//planar().
bool MGCylinder::flat(
	const MGBox& uvbox,
	double tol,		//Tolerance allowed to regard flat
					//(Allowed distance from a Cylinder).
	int& direction,	//   1: u-direction is more non flat.
					//   0: v-direction is more non flat.
	MGPosition& P,	//Position of the flat plane will be output.
	MGUnit_vector& N//Normal of the flat plane will be output.
)const{
	direction=1;
	N=normal(uvbox.mid());
	const MGInterval& urange=uvbox[0];
	const MGInterval& vrange=uvbox[1];
	MGEllipse el(mgORIGIN_2D,
		MGVector(m_ellipse.major_len(),0.),
		MGVector(0.,m_ellipse.minor_len()),
		urange);
	const MGBox& elbx=el.box();
	P=elbx.mid();
	P+=m_ellipse.center();
	P+=m_axis.eval(vrange.mid_point());
	double a=elbx[0].length().value(), b=elbx[1].length().value();
	if(a<b) a=b;
	double uspan=el.m_prange[1]-el.m_prange[0];
	if(uspan<mgPAI){
		MGMatrix mat;mat.set_x_axis(el.end_point()-el.start_point());
		el*=mat;
		const MGBox& elbx2=el.box();
		double a2=elbx2[0].length().value(), b2=elbx2[1].length().value();
		if(a2<b2) a2=b2;
		if(a>a2) a=a2;
	}

	return a<=tol;
}

//This is the same as flat except that this does not have the arguments P, N.
//***** the flatness is tested only approximately. This is for exclusive use of
//tessellation.
bool MGCylinder::flat_tess(
	double u0,double u1,	//u range from u0 to u1.
	double v0,double v1,	//v range from v0 to v1.
	double tol,		//Tolerance allowed to regart flat
					//(Allowed distance from a Cylinder).
	bool& direction,//   1: u-direction is more non flat.
					//   0: v-direction is more non flat.
	double max_edge_len
)const{
	direction=true;
	double um=(u0+u1)*0.5;

	size_t i;	//id of Pn[].
	MGVector Pn[3], Nn[3];
	Pn[0]=eval(u0, v0); Nn[0]=normal(u0,v0);
	Pn[1]=eval(um, v0); Nn[1]=normal(um,v0);
	Pn[2]=eval(u1, v0); Nn[2]=normal(u1,v0);

	MGVector P=Pn[0]; MGVector VN=Nn[0]; 
	for(i=1; i<3; i++){P+=Pn[i]; VN+=Nn[i];}
	P/=3.;
	MGUnit_vector N(VN);

	double x, d=P%N;
	double dist[3];
	bool is_flat=true;
	for(i=0; i<3; i++){
		x=dist[i]=d-Pn[i]%N;
		if(x<0.) x=-x;
		if(x>tol) is_flat=false;
	}
	if(is_flat){
		MGVector P02(Pn[2],Pn[0]);
		double lenu=P02%P02;
		MGVector axisV=m_axis.eval(v0)-m_axis.eval(v1);
		double lenv=axisV%axisV;
		direction=lenu>=lenv;
		if(max_edge_len<=0.) return true;
		double melen2=max_edge_len;
		melen2*=melen2;
		if(direction){
			if(lenu>melen2) return false;
		}else{
			if(lenv>melen2) return false;
		}
		return true;
	}
	return false;
}

bool MGCylinder::in_range(double u, double v) const{
	return m_ellipse.in_range(u) && m_axis.in_range(v);
}

// Surface と Curve の交点を求める。
//Compute intesection of Cylinder and Curve.
MGCSisect_list MGCylinder::isect(const MGCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Cylinder and Curve.
MGCSisect_list MGCylinder::isect(const MGRLBRep& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Cylinder and Curve.
MGCSisect_list MGCylinder::isect(const MGEllipse& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Cylinder and Curve.
MGCSisect_list MGCylinder::isect(const MGLBRep& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Cylinder and Curve.
MGCSisect_list MGCylinder::isect(const MGSurfCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Curve の交点を求める。
//Compute intesection of Cylinder and Curve.
MGCSisect_list MGCylinder::isect(const MGBSumCurve& curve) const{
	return curve.isect(*this);
}

// Surface と Surface の交線を求める。
//Compute intesection of Sphere and Surface.
MGSSisect_list MGCylinder::isect(const MGSurface& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGCylinder::isect(const MGPlane& srf2) const{
	return intersectPl(srf2);
}
MGSSisect_list MGCylinder::isect(const MGSphere& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGCylinder::isect(const MGCylinder& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGCylinder::isect(const MGSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGCylinder::isect(const MGRSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGCylinder::isect(const MGBSumSurf& srf2) const{
	return intersect(srf2);
}

//isect_direction() is used by isect_startPt() to define which constant
//parameter line should be used to compute intersection, and what
//incremental value be used for the parameter.
//Function's return value is direction to get next intersection(with dt).
//When =1: u=const direction, =0: v=const, =-1: cannot get intersection.
int MGCylinder::isect_direction(
	const MGFSurface& sf2,	//Second surface for the intersection.
	size_t m1,		//id of uvuvS that indicates this surface's parameter
		//position in uvuvS. (uvuvS(m1), uvuvS(m1+1))=(u,v) of this surface.
	MGPosition& uvuvS,//start parameter (u,v) pair of this surface and sf2.
	double& du,	//Incremental value of the parameter kind of kdt will be output.
	double& dv, //Right dt will be output according to the function's output =0,1.
	double acuRatio		//acuracy ratio.
)const{
	const MGPlane* pl=dynamic_cast<const MGPlane*>(&sf2);
	if(!pl)
		return MGFSurface::isect_direction(sf2,m1,uvuvS,du,dv,acuRatio);

	double du2,dv2;
	isect_dt(uvuvS[m1],uvuvS[m1+1],du2,dv2,acuRatio);
	du=du2;
	dv=dv2;
	const MGUnit_vector& pln=pl->normal();
	MGUnit_vector axisv=m_axis.direction();
	double angl=pln%axisv;
	if(fabs(angl)<.01) return 0;
	return 1;
}

#define MGCylinder_isect_div_num1 12.
#define MGCylinder_isect_div_num2 5.
//isect_dt computes incremental values of u and v direction for the intersection
//computation at parameter position (u,v).
void MGCylinder::isect_dt(
	double u, double v, double& du, double& dv,
	double acuRatio		//acuracy ratio.
)const{
	double alfa=isect_dt_coef(0);
	du=m_ellipse.param_span()/MGCylinder_isect_div_num1*alfa*acuRatio;
	dv=m_axis.param_span()/MGCylinder_isect_div_num2*alfa*acuRatio;
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGCylinder::isect_incr_pline(
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
		u=m_ellipse.range(uv[0]+du); v=uv[1];
		//Compute u=const v-parameter line of this surface in pline.
		return parameter_curve(kdt,u);
	}else{
		v=m_axis.range(uv[1]+dv); u=uv[0];
		//Compute v=const u-parameter line in pline.
		return parameter_curve(kdt,v);
	}
}

//"isect_inner_dt" is a dedicated function of isect_startPt,
// comutes adequate incremental parameter value(du,dv) and parameter line kind
//kdt(u=const or v=const).
void MGCylinder::isect_inner_dt(
	size_t n,	//num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,//intersection point obtained last(of this).
	double& du, double& dv,	//incremental length from previous to uvnow is input.
				//New du or dv will be output according to kdt's return value.
	int& kdt,	//Parameter kind used so far is input, will be output as:
				//=1:parameter line kind(u=const), =0: v=const,
				//=-1:should halt computation since incremental value is zero.
	double acuRatio	//Accurate ratio.
) const{
	double uerr=param_error_u();
	double verr=param_error_v();
	double abdu=fabs(du), abdv=fabs(dv);
	double ratio=acuRatio;
	if(abdu<uerr){
		if(abdv<verr){kdt=-1; return;}
		else kdt=0;
	}else if(abdv<verr) kdt=1;
	else{
		MGVector dfdu=eval(uvnow,1,0), dfdv=eval(uvnow,0,1);
		MGVector dfdt=dfdu*du+dfdv*dv;
		double fuu=dfdu%dfdu, fuv=dfdu%dfdv;
		double dffu=fuu*du+fuv*dv;
		double cos2=dffu*dffu/(dfdt%dfdt)/(dfdu%dfdu);
		double rcz=MGTolerance::rc_zero();
		rcz*=5.;
		rcz*=rcz;
		if(cos2<=rcz) kdt=0; else kdt=1;
	}

//Define new dt, kdt.
	double dt,dtold;
	if(kdt){
		dtold=du;
		dt=m_ellipse.param_span()/MGCylinder_isect_div_num1;
	}else{
		dtold=dv;
		dt=m_axis.param_span()/MGCylinder_isect_div_num2;
	}
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

typedef MGPosition_list::iterator plitr;
bool uvcompare(plitr i1,plitr i2){return (*i1)[0]<(*i2)[0];};

//Compute the intersection line of this and the plane pl.
MGSSisect_list MGCylinder::intersectPl(const MGPlane& pl)const{
	MGSSisect_list lst(this,&pl);
	MGPosition_list uvuv_list2;
	double error=MGTolerance::set_wc_zero(MGTolerance::line_zero()*.5);
		//Save the error.
	isect_boundary(pl,uvuv_list2);
	MGTolerance::set_wc_zero(error);//Restore the error.
	size_t numi=uvuv_list2.size();
	if(numi){
	
	size_t j;
	std::vector<plitr> iss(numi);
	plitr i=uvuv_list2.begin(), ie=uvuv_list2.end();
	for(j=0; j<numi; j++,i++) iss[j]=i;
	std::sort(iss.begin(), iss.end(), uvcompare);
	MGPosition_list uvuv_list;
	for(j=0; j<numi; j++,i++) uvuv_list.append(*(iss[j]));
	//cout<<uvuv_list;
	double uspan=m_ellipse.param_span();
	double verr=m_axis.param_error();
	while(uvuv_list.entries()>1){
		plitr pi=uvuv_list.begin();
		MGPosition& uvuv1=*(pi++);
		MGPosition& uvuv2=*pi;
		if(MGREqual_base(uvuv1[0], uvuv2[0], uspan)){
			MGCurve* iline=
				new MGStraight(pl.eval(uvuv2[2],uvuv2[3]),pl.eval(uvuv1[2],uvuv1[3]));
			double t1=iline->param_s(), t2=iline->param_e();
			MGCurve* param1=
				new MGStraight(MGPosition(uvuv2[0],uvuv2[1]), MGPosition(uvuv1[0],uvuv1[1]));
			param1->change_range(t1,t2);
			MGCurve* param2=
				new MGStraight(MGPosition(uvuv2[2],uvuv2[3]), MGPosition(uvuv1[2],uvuv1[3]));
			param2->change_range(t1,t2);
			lst.append(iline,param1,param2);
			uvuv_list.m_Plist.pop_front();
			uvuv_list.m_Plist.pop_front();
		}else if(fabs(uvuv1[1]-uvuv2[1])<=verr){
			MGCurve* iline=
				new MGEllipse(m_ellipse+m_axis.eval(uvuv1[1]));
			double t1=iline->param_s(), t2=iline->param_e();
			MGCurve* param1=
				new MGStraight(MGPosition(uvuv2[0],uvuv2[1]), MGPosition(uvuv1[0],uvuv1[1]));
			param1->change_range(t1,t2);
			MGPosition elcenter=iline->eval((t1+t2)*.5);
			MGPosition pluvc;
			pl.on(elcenter,pluvc);
			MGCurve* param2=
				new MGEllipse(
				MGPosition(uvuv1[2],uvuv1[3]),	//Starting point
				pluvc,
				MGPosition(uvuv2[2],uvuv2[3]));//Ending point
			param2->change_range(t1,t2);
			lst.append(iline,param1,param2);
			uvuv_list.m_Plist.pop_front();
			uvuv_list.m_Plist.pop_front();
		}else{
			lst.append(isect_with_surf(uvuv_list,pl));
		}
	}

	}
	return lst;
}

//Intersection of Surface and a straight line.
MGCSisect_list MGCylinder::isectSl(
	const MGStraight& sl,
	const MGBox& uvbox //indicates if this surface is restrictied to the parameter
					//range of uvbox. If uvbox.is_null(), no restriction.
)const{
	assert(m_ortho);//currently, m_axis must be orthogonal to m_ellipse.
	//cout<<uvbox<<*this<<endl;

	const MGVector& axisdir=m_axis.direction();
	MGMatrix mat;mat.set_axis(axisdir,2);
	MGEllipse el2(2,m_ellipse*mat);//cout<<el2;
	MGStraight sl21(2,sl*mat);//cout<<sl21<<endl;
	MGCCisect_list isects=el2.isect(sl21);//cout<<isects<<endl;

	MGCSisect_list list(&sl,this);
	if(!isects.entries()) return list;

	const MGVector& sldir=sl.direction();
	MGVector slbyaxis=sldir*axisdir;
	MGMatrix mat2; mat2.set_axis(slbyaxis,2);
	MGStraight sl22(2,sl*mat2);//cout<<sl22;
	MGCCisect_list::CCiterator i=isects.begin(), ie=isects.end();
	for(; i!=ie; i++){
		double u=(*i).param1(), t1=(*i).param2();
		MGStraight* slcyl=static_cast<MGStraight*>(parameter_curve(true,u));
		MGStraight slcyl2(2,(*slcyl)*mat2);//cout<<slcyl2;
		delete slcyl;
		MGCCisect visect; slcyl2.relation(sl22,visect);//cout<<visect;
		double v=visect.param1(), t2=visect.param2();
		double t=(t1+t2)*.5;//To increase the accuracy.
		if(!uvbox.is_null()){
			if(!uvbox[0].includes(u)) continue;
			if(!uvbox[1].includes(v)) continue;
		}
		list.append(sl.eval(t),t,MGPosition(u,v));
	}
	return list;
}

// 柱面を反転する。ノーマルを逆方向にする。
//Negate the normal of the Cylinder.
void MGCylinder::negate(
		int is_u)// Negate along u-direction if is_u is ture,
				// else along v-direction.
{
	if(is_u) m_ellipse.negate();
	else m_axis.negate();
}

//Obtain parameter value if this surface is negated by "negate()".
//Negate along u-direction if is_u is ture,
// else along v-direction.
MGPosition MGCylinder::negate_param(const MGPosition& uv, int is_u)const{
	MGPosition uvnew(uv);
	if(is_u) uvnew(0)=m_ellipse.negate_param(uv[0]);
	else     uvnew(1)=m_axis.negate_param(uv[1]);
	return uvnew;
}

//C1連続曲面の一定オフセット関数
//オフセット方向は、ノーマル方向を正とする。トレランスはline_zero()を使用している。
//戻り値は、オフセット面のオートポインターが返却される。
//Surface offset. positive offset value is offset normal direction.
//the radius of curvature must be larger than offset value.
//line_zero() is used.
std::auto_ptr<MGSurface> MGCylinder::offset_c1(
	double ofs_value,	//オフセット量
	int& error//エラーコード 0:成功 -1:面におれがある*
			// -2:曲率半径以上のオフセット不可 -3:面生成コンストラクタエラー
)const{
	MGVector S=m_ellipse.eval(m_ellipse.param_s(),1);
	const MGVector& M=m_ellipse.major_axis();
	const MGVector& A=m_axis.direction();
	if((S*A)%M <0.) ofs_value*=-1.;
	double r1=m_ellipse.m_r;
	double r2=r1+ofs_value;
	if(r2<=0.){
		error=-2;
		return std::auto_ptr<MGSurface>();
	}
	MGEllipse el2(m_ellipse);
	double scale=r2/r1;
	el2.m_m*=scale;
	el2.m_n*=scale;
	el2.m_r*=scale;
	el2.update_mark();
	std::auto_ptr<MGSurface> surf(new MGCylinder(el2,m_axis));
	return surf;
}

// 指定点が面上にあるか調べる。（面上ならばtrue）
//Test if a point is on the Cylinder. If on the Cylinder, return true.
bool MGCylinder::on(
	const MGPosition& point,	//A point to test 指定点
	MGPosition& puv				//Parameter value of the Cylinder will be
								//returned.
)const{
	puv=closest(point);
	return ((point-eval(puv)).len()<=MGTolerance::wc_zero());
}

//Test if input (u,v) is parameter value on a perimeter of the surface.
//If u or v is on a perimeter, (u,v) will be updated to the perimeter value.
bool MGCylinder::on_a_perimeter(
	double& u, double& v,	//Surface parameter (u,v)
	size_t& perim_num	//if function returns true,
						//the perimete rnumber is output.
)const{
	if(!in_range(u,v)) return false;
	MGSTRAIGHT_TYPE slt=m_axis.straight_type();
	if(slt!=MGSTRAIGHT_UNLIMIT){
		double vspan=m_axis.param_span();
		double v0=m_axis.param_s();
		if(MGREqual_base(v,v0,vspan)){
			perim_num=0; v=v0;
			return true;
		}
		if(slt!=MGSTRAIGHT_HALF_LIMIT){
			double v1=m_axis.param_e();
			if(MGREqual_base(v,v1,vspan)){
				perim_num=2; v=v1;
				return true;
			}
		}
	}
	double uspan=m_ellipse.param_span();
	double u0=m_ellipse.param_s();
	if(MGREqual_base(u,u0,uspan)){
		perim_num=3; u=u0;
		return true;
	}
	double u1=m_ellipse.param_e();
	if(MGREqual_base(u,u1,uspan)){
		perim_num=1; u=u1;
		return true;
	}
	return false;
}

//Obtain parameter space error.
double MGCylinder::param_error() const{
	double uerror=m_ellipse.param_error();
	double verror=m_axis.param_error();
	return sqrt(uerror*uerror+verror*verror);
}

// パラメータ範囲を返す。
//Return parameter range of the Cylinder(Infinite box).
MGBox MGCylinder::param_range() const{
	return MGBox(m_ellipse.param_range(), m_axis.param_range());
}

// Compute parameter curve.
//Returned is newed area pointer, and must be freed by delete.
MGCurve* MGCylinder::parameter_curve(
	int is_u				//Indicates x is u-value if is_u is true.
	, double x				//Parameter value.
							//The value is u or v according to is_u.
)const{
	if(is_u){
		MGStraight* sl=new MGStraight(m_axis);
		sl->m_root_point=m_ellipse.eval(x);
		return sl;
	}else{
		MGEllipse* el=new MGEllipse(m_ellipse);
		el->m_center+=m_axis.eval(x);
		return el;
	}
}

//Compute part of the surface limitted by the parameter range bx.
//bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
MGCylinder* MGCylinder::part(
	const MGBox& uvbx,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
)const{
	MGCylinder* cyl=new MGCylinder(*this);
	cyl->m_ellipse.limit(uvbx[0]);
	cyl->m_axis.limit(uvbx[1]);
	return cyl;
}

// i must be < perimeter_num().
//When perimeter_num()==0, this function is undefined.
MGCurve* MGCylinder::perimeter_curve(size_t i) const{
	if(i==0){
		return new MGEllipse(m_ellipse+m_axis.start_point());
	}else if(i==2){
		return new MGEllipse(m_ellipse+m_axis.end_point());
	}else if(i==1){
		return new MGStraight(m_axis+m_ellipse.end_point());
	}else{
		return new MGStraight(m_axis+m_ellipse.start_point());
	}
}

// 与えられた点にもっとも近い面上の垂直のパラメータ値を返却する。
//Return the nearest perpendicular point of the Cylinder from P.
// Function's return value is whether point is obtained(1) or not(0)
int MGCylinder::perp_point(
	const MGPosition& P,// 与えられた点
	MGPosition& uv,		//Parameter value of the Cylinder will be output
	const MGPosition* uvguess	// guess
)const{
	assert(m_ortho);
	double v;
	int obtained;
	int range=vrange(P,v);
	MGEllipse el(m_ellipse);
	if(range==1){
		el+=m_axis.end_point();
		obtained=0;
	}else if(range==-1){
		el+=m_axis.start_point();
		obtained=0;
	}else{
		el+=m_axis.eval(v);
		obtained=1;
	}
	double u=el.closest(P);
	uv=MGPosition(u,v);
	return obtained;
}

//Return all(actually one) foots of perpendicular straight lines from P.
MGPosition_list MGCylinder::perps(
	const MGPosition& P				// Point of a space(指定点)
) const{
	assert(m_ortho);
	MGPosition_list list;
	double v;
	int vinrange=vrange(P,v);
	if(vinrange==0){
		MGPosition p2d=
			P*MGTransf(m_ellipse.major_axis(), m_ellipse.minor_axis(),axis_point(v));
		double theta[4];		// Intersection points are maixmum 4.
		size_t nump=m_ellipse.perp2d(p2d.ref(0), p2d.ref(1), theta);
		for(size_t i=0; i<nump; i++){
			double u=m_ellipse.radian_to_gp(theta[i]);
			if(m_ellipse.in_range(u)) list.append(MGPosition(u,v));
		}
		return list;
	}else return list;
}

// 入力パラメータをパラメータ範囲でまるめて返却する。
//Round the input uv into parameter range of the Cylinder, 
//return the same value as input.
MGPosition MGCylinder::range(const MGPosition& uv) const{
	return MGPosition(m_ellipse.range(uv[0]), m_axis.range(uv[1]));
}

//Return the space dimension.
size_t MGCylinder::sdim() const{return 3;}

//Obtain the v parameter value of the neareast point from the point P to the axis.
//Function's return value is if the parameter value v is in the range of this
//cylinder: -1: below the range, 0:in the range, 1:above the range.
int MGCylinder::vrange(const MGPosition& P, double& v)const{
	MGStraight sl(m_axis);
	sl+=m_ellipse.center();
	v=sl.perp_param(P);
	MGInterval vrng=sl.param_range();
	int obtained;
	if(vrng<v){
		v=vrng[1].value();
		obtained=1;
	}else if(vrng>v){
		v=vrng[0].value();
		obtained=-1;
	}else{
		obtained=0;
	}
	return obtained;
}

//メンバデータを読み込む関数
void MGCylinder::ReadMembers(MGIfstream& buf){
	MGSurface::ReadMembers(buf);
	m_ellipse.ReadMembers(buf);
	m_axis.ReadMembers(buf);
}

//メンバデータを書き込む関数
void MGCylinder::WriteMembers(MGOfstream& buf) const{
	MGSurface::WriteMembers(buf);
	m_ellipse.WriteMembers(buf);
	m_axis.WriteMembers(buf);
}
