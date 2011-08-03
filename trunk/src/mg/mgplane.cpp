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
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/Plane.h"
#include "mg/SPointSeq.h"
#include "mg/SBRep.h"
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

// コンストラクタ
// 初期化なしで平面を生成する。
MGPlane::MGPlane():MGSurface(),m_uknotV(0),m_vknotV(0){;}

//Copy constructor.
MGPlane::MGPlane(const MGPlane& pl)
:MGSurface(pl),m_d(pl.m_d),m_normal(pl.m_normal),
m_uderiv(pl.m_uderiv),m_vderiv(pl.m_vderiv),
m_root_point(pl.m_root_point),m_uknotV(0),m_vknotV(0){;}

// Construct a plane by changing this space dimension or
// ordering the coordinates.
MGPlane::MGPlane(
	size_t dim,				// New space dimension.
	const MGPlane& plane,	// Original Surface B-rep.
	size_t start1,	 		// Destination order of new Surface.
	size_t start2)			// Source order of original Surface.
:MGSurface(plane),
	m_uderiv(dim,plane.m_uderiv,start1,start2),
	m_vderiv(dim,plane.m_vderiv,start1,start2),
	m_root_point(dim,plane.m_root_point,start1,start2)
	,m_uknotV(0),m_vknotV(0){
	update_mark();
	MGUnit_vector U=m_uderiv,V;
	U.orthonormal(m_vderiv, V, m_normal);
	m_d=m_normal%m_root_point;
}

// Construct a plane from the coefficients of the pla equation:
//     a*x+b*y+c*z=d.
//     coefficients a,b,c,d are provided by double array g[4].
MGPlane::MGPlane(
	const double g[4],		//coefficients g[0]=a, b=g[1], c=g[2], d=g[3].
	const double* root_point//When root_point!=0, root_point[.] are (x,y,z) values of the root point.
):MGSurface(),m_uknotV(0),m_vknotV(0){
	MGVector v(3,g);
	m_normal=v;
	double vlen=v.len();
	if(MGMZero(vlen))
		m_d=0.;
	else
		m_d=g[3]/vlen;
	m_root_point=m_normal*m_d;
	MGUnit_vector U,V;
	m_normal.orthonormal(m_normal, U, V);
	m_uderiv=U; m_vderiv=V;
	if(root_point){//When root point is specified.
		m_root_point=eval(uv(MGPosition(3,root_point)));
	}
}

// 平面のノーマルと平面の原点からの距離を指定して面を生成する。
MGPlane::MGPlane(
	const MGUnit_vector& normal,	//Normal of the plane.
	double d				//distance from origin of the plane.
							//When normal=(a,b,c,...), d=a*x+b*y+c*z+.... .
):MGSurface(), m_normal(normal), m_d(d), m_root_point(normal*d)
	,m_uknotV(0),m_vknotV(0){
	MGUnit_vector U,V;
	m_normal.orthonormal(m_normal, U, V);
	m_uderiv=U; m_vderiv=V;
}

// 点と平面のノーマルを指定して面を生成する。
MGPlane::MGPlane (
	const MGUnit_vector& norml,
	const MGPosition& p
):MGSurface(), m_normal(norml), m_d(p%norml), m_root_point(p)
	,m_uknotV(0),m_vknotV(0){
	MGUnit_vector U,V;
	m_normal.orthonormal(m_normal, U, V);
	m_uderiv=U; m_vderiv=V;
}

// 直線と直線に乗らない点を指定して面を生成する。
MGPlane::MGPlane(
	const MGStraight& st,
	const MGPosition& point
):MGSurface(), m_root_point(point), m_uderiv(st.direction().normalize())
	,m_uknotV(0),m_vknotV(0){
	MGUnit_vector U,V;
	U=m_uderiv;
	U.orthonormal(st.root_point()-m_root_point, V, m_normal);
	m_vderiv=V;
	m_d=m_normal%m_root_point;
}

// 点、ｕ方向、ｖ方向を指定して平面を生成する。
MGPlane::MGPlane (
	const MGVector& u,
	const MGVector& v,
	const MGPosition& p
):MGSurface(), m_root_point(p) ,m_uderiv(u),m_vderiv(v)
	,m_uknotV(0),m_vknotV(0){
	MGUnit_vector uunit=u,vunit, normal;
	uunit.orthonormal(v, vunit, normal);
	m_normal=normal;
	m_d=p%m_normal;
}

// Construct a plane by interpolating two planes.
// If two planes intersect, interpolation is rotation.
//If two planes are parallel, interpolation is parallel move.
MGPlane::MGPlane(
	const MGPlane& plane1,
	const MGPlane& plane2,
	double t				//Input ratio.
				// When t=0, the plane will be plane1.
				// When t=1, the plane will be plane2.
):MGSurface(plane1),m_uknotV(0),m_vknotV(0){
	update_mark();
	double onemt=1.-t;
	if((plane1.normal()).parallel(plane2.normal())){
		m_normal=plane1.normal();
		m_d=onemt*plane1.distance()+t*plane2.distance();
		m_uderiv=plane1.u_deriv(); m_vderiv=plane1.v_deriv();
	}
	else{
		const MGVector& e1=plane1.normal();
		const MGVector& e2=plane2.normal();
		double ratio[2];
		m_normal=e1.interpolate_by_rotate(t,e2,ratio);
	    m_d=plane1.distance()*ratio[0] + plane2.distance()*ratio[1];

		MGUnit_vector uunit,vunit;
		m_normal.orthonormal(m_normal, uunit, vunit);
		m_uderiv=uunit; m_vderiv=vunit;
	}
	m_root_point=m_normal*m_d;
}

//Construct a plane from three points on the plane.
MGPlane::MGPlane(
	const MGPosition& P1,
	const MGPosition& P2,
	const MGPosition& P3
): MGSurface(), m_uknotV(0),m_vknotV(0){
	MGPosition P((P1+P2+P3)/3.);
	MGVector v1(P1,P);
	MGVector v2(P2,P);
	MGVector v3(P3,P);
	double a1=v1%v2; if(a1<0.) a1=-a1;
	double a2=v2%v3; if(a2<0.) a2=-a2;
	double a3=v3%v1; if(a3<0.) a3=-a3;
	MGVector N;
	if(a1<a2)
		if(a1<a3) N=v1*v2;
		else N=v3*v1;
	else if(a3<a2) N=v3*v1;
	else N=v2*v3;
	*this=MGPlane(N, P);
}

//////////Destructor//////////////
MGPlane::~MGPlane(){
	if(m_uknotV) delete m_uknotV;
	if(m_vknotV) delete m_vknotV;
}

// メンバ関数

//Gets parameters(a,b,c,d) of the plane expression a*x+b*y+c*z=d.
void MGPlane::abcd(double g[4]) const{
//g[.]=(a,b,c,d)
	for(size_t i=0; i<3; i++) g[i]=m_normal.ref(i);
	g[3]=m_d;
}

// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox MGPlane::box_limitted(
	const MGBox& uv	// Parameter Range of the surface.
) const{
	const MGInterval &urng=uv[0], &vrng=uv[1];
	if(urng.finite() && vrng.finite()){
		double u0=urng.low_point(), u1=urng.high_point();
		double v0=vrng.low_point(), v1=vrng.high_point();
		MGBox bx(eval(u0,v0),eval(u1,v1));
		bx.expand(eval(u0,v1));
		bx.expand(eval(u1,v0));
		return bx;
	}

	double um,vm;
	if(urng.finite()){
		um=urng.mid_point();
		if(vrng.finite_above()) vm=vrng.high_point();
		else if(vrng.finite_below()) vm=vrng.low_point();
		else vm=0.;
	}else{
		vm=vrng.mid_point();
		if(urng.finite_above()) um=urng.high_point();
		else if(urng.finite_below()) um=urng.low_point();
		else um=0.;
	}

	MGPosition origin(eval(um,vm));
	MGStraight uline;
	uline.set_straight(MGSTRAIGHT_UNLIMIT, m_uderiv, origin);
	MGStraight vline;
	vline.set_straight(MGSTRAIGHT_UNLIMIT, m_vderiv, origin);
	uline.limit(uv.ref(0));
	vline.limit(uv.ref(1));

	MGSTRAIGHT_TYPE ut=uline.straight_type();
	if(ut==MGSTRAIGHT_EMPTY) return vline.box();
	MGSTRAIGHT_TYPE vt=vline.straight_type();
	if(vt==MGSTRAIGHT_EMPTY) return uline.box();

	if(ut==MGSTRAIGHT_SEGMENT){
		MGVector uvec=uline.end_point()-uline.start_point();
		return ((vline+uvec).box() | vline.box());
	}else if(vt==MGSTRAIGHT_SEGMENT){
		MGVector vvec=vline.end_point()-vline.start_point();
		return ((uline+vvec).box() | uline.box());
	}else return uline.box() | vline.box();
}

//Changing this object's space dimension.
MGPlane& MGPlane::change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	m_uderiv=MGVector(sdim,m_uderiv,start1,start2);
	m_vderiv=MGVector(sdim,m_vderiv,start1,start2);
	m_root_point=MGPosition(sdim,m_root_point,start1,start2);
	MGUnit_vector U=m_uderiv,V;
	U.orthonormal(m_vderiv, V, m_normal);
	m_d=m_normal%m_root_point;
	update_mark();
	return *this;
}

//Return minimum box that includes whole of the surface.
MGBox* MGPlane::compute_box() const{
	MGInterval span(MGINTERVAL_INFINITE);
	MGBox bx(span,span); //u span, and v span;
	MGBox plbox=box_limitted(bx);
	MGBox* bxp=new MGBox(plbox);
	return bxp;
}

//Compute the closest point parameter value (u,v) of this surface
//from a point.
MGPosition MGPlane::closest(const MGPosition& point) const{
	return param(point);
}

//Construct new surface object by copying to newed area.
//User must delete this copied object by "delete".
MGPlane* MGPlane::clone() const{return new MGPlane(*this);}

//Construct new surface object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGPlane* MGPlane::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2 		// Source order of this line.
)const{
	return new MGPlane(sdim,*this,start1,start2);
}

// 自身と与えられた点との距離を返却する。
double MGPlane::distance(const MGPosition& p) const{
	return  m_d - m_normal%p;
}

//Evaluate surface data.
MGVector MGPlane::eval(
	double u, double v	// Parameter value of the surface.
	, size_t ndu			// Order of derivative along u.
	, size_t ndv			// Order of derivative along v.
	) const{
	MGVector data;
	if(ndu==0){
		if(ndv==0)
			data=m_root_point+m_uderiv*u+m_vderiv*v;
		else if(ndv==1) data=m_vderiv;
	}
	else if(ndu==1 && ndv==0) data=m_uderiv;
	else data=MGVector(0.,0.,0.);

	return data;
}

//Evaluate right continuous surface data.
//Evaluate all positional data, 1st and 2nd derivatives.
void MGPlane::eval_all(
	double u, double v,		// Parameter value of the surface.
	MGPosition& f,			// Positional data.
	MGVector&   fu,			// df(u,v)/du
	MGVector&   fv,			// df/dv
	MGVector&   fuv,		// d**2f/(du*dv)
	MGVector&   fuu,		// d**2f/(du**2)
	MGVector&   fvv			// d**2f/(dv**2)
	) const{
	f=m_root_point+m_uderiv*u+m_vderiv*v;
	fu=m_uderiv; fv=m_vderiv;
	fuv=fuu=fvv=MGVector(0.,0.,0.);
}

// Exchange parameter u and v.
MGSurface& MGPlane::exchange_uv(){
	m_normal = -m_normal;
	MGVector save(m_uderiv);
	m_uderiv=m_vderiv;
	m_vderiv=save;
	return *this;
}

//Test if the surface is flat or not within the parameter value rectangle of uvbox.
//Function's return value is:
//	true: if the surface is flat
//  false: if the surface is not falt.
//When this is not falt, the direction that indicates which direction the surface
//should be divided will be output.
//***** the flatness is tested only approximately. This is for exclusive use of
//planar().
bool MGPlane::flat(
	const MGBox& uvbox,
	double tol,		//Tolerance allowed to regard flat
					//(Allowed distance from a plane).
	int& direction,	//   1: u-direction is more non flat.
					//   0: v-direction is more non flat.
	MGPosition& P,	//Position of the flat plane will be output.
	MGUnit_vector& N//Normal of the flat plane will be output.
)const{
	direction=1;
	P=m_root_point;
	N=m_normal;
	return true;
}

//This is the same as flat except that this does not have the arguments P, N.
//***** the flatness is tested only approximately. This is for exclusive use of
//tessellation.
bool MGPlane::flat_tess(
	double u0,double u1,	//u range from u0 to u1.
	double v0,double v1,	//v range from v0 to v1.
	double tol,		//Tolerance allowed to regart flat
					//(Allowed distance from a plane).
	bool& direction	//   1: u-direction is more non flat.
					//   0: v-direction is more non flat.
)const{
	MGVector P00=eval(u0,v0), P01=eval(u0,v1),
			P10=eval(u1,v0), P11=eval(u1,v1);
	MGVector alongU0=P10-P00;
	MGVector alongU1=P11-P01;
	double lenu=alongU0%alongU0+alongU1%alongU1;
	MGVector alongV0=P01-P00;
	MGVector alongV1=P11-P10;
	double lenv=alongV0%alongV0+alongV1%alongV1;
	direction=lenu>=lenv;
	return true;
}

//Return virtual knot value of plane.
double MGPlane::knot_u(size_t i)const {
	assert(i<=1);
	if(i) return param_e_u(); else return param_s_u();
}

double MGPlane::knot_v(size_t i)const {
	assert(i<=1);
	if(i) return param_e_v(); else return param_s_v();
}

//Returns the u knot vector.
const MGKnotVector& MGPlane::knot_vector_u() const{
	if(!m_uknotV){
		m_uknotV=new MGKnotVector(1,1);
		(*m_uknotV)(0)=param_s_u(); (*m_uknotV)(1)=param_e_u();
	}
	return *m_uknotV;
}
MGKnotVector& MGPlane::knot_vector_u(){
	const MGPlane* cpl=this;
	const MGKnotVector& ckt=cpl->knot_vector_u();
	MGKnotVector* kt=const_cast<MGKnotVector*>(&ckt);
	return *kt;
}

//Returns the v knot vector.
const MGKnotVector& MGPlane::knot_vector_v() const{
	if(!m_vknotV){
		m_vknotV=new MGKnotVector(1,1);
		(*m_vknotV)(0)=param_s_v(); (*m_vknotV)(1)=param_e_v();
	}
	return *m_vknotV;
}
MGKnotVector& MGPlane::knot_vector_v(){
	const MGPlane* cpl=this;
	const MGKnotVector& ckt=cpl->knot_vector_v();
	MGKnotVector* kt=const_cast<MGKnotVector*>(&ckt);
	return *kt;
}

// 平面を反転する。ノーマルを逆方向にする。
void MGPlane::negate(
			int is_u)	// Negate along u-direction if is_u is ture,
						// else along v-direction.
{
	m_normal = -m_normal;
	if(is_u) m_uderiv = -m_uderiv;
	else     m_vderiv = -m_vderiv;
	m_d = -m_d;
}

//Obtain parameter value if this surface is negated by "negate()".
// Negate along u-direction if is_u is ture,
// else along v-direction.
MGPosition MGPlane::negate_param(const MGPosition& uv, int is_u)
const{
	double u=uv(0), v=uv(1);
	if(is_u) u=-u; else v=-v;
	return MGPosition(u,v);
}

// 指定点が面上にあるか調べる。（面上ならばtrue）
//Test if a point is on the plane. If on the plane, return true.
bool MGPlane::on(const MGPosition& point)const{	
	return MGAZero(distance(point));
}

// 与えられた誤差内で点が面上にあるかどうかテストする。
bool MGPlane::on (
	const MGPosition& point,		// 指定点                        
	MGPosition& puv				// Parameter value will be returned.                    
)const{	
	double t=distance(point);
	puv=uv(point);
	return MGAZero(t);
}

// 直線が平面上にあるか調べる。（平面上ならばtrue）
bool MGPlane::on(
	const MGStraight& sl	// 直線
	) const{
	MGPosition uv(2);
	return sl.direction().orthogonal(normal()) && on(sl.root_point(), uv);
}

//Return plane parameter value of a point on the plane. 
//If input point is not on the plane, returned is
//the nearest point parameter of the plane.
MGPosition MGPlane::param(
	const MGPosition& p		// 指定点
	) const
{	return uv(p);}

//Obtain parameter space tolerance==rc_zero();
double MGPlane::param_error() const{
	double ue=param_error_u(), ve=param_error_v();
	return sqrt(ue*ue+ve*ve);
}
double MGPlane::param_error_u() const{return m_uderiv.len()*MGTolerance::wc_zero();}
double MGPlane::param_error_v() const{return m_vderiv.len()*MGTolerance::wc_zero();}

// パラメータ範囲を返す。
MGBox MGPlane::param_range() const{
	MGInterval iv(MGINTERVAL_INFINITE);
	MGBox box(iv,iv);
	return box;
}

// Compute parameter curve.
//Returned is newed area pointer, and must be freed by delete.
MGCurve* MGPlane::parameter_curve(
	int is_u				//Indicates x is u-value if is_u is true.
	, double x				//Parameter value.
							//The value is u or v according to is_u.
	) const
{
	MGPosition origin;
	const MGVector* dir; const MGVector* dir2;

	if(is_u){
		dir2=&u_deriv();
		dir=&v_deriv();
	}else{
		dir2=&v_deriv();
		dir=&u_deriv();
	}
	origin=root_point()+x*(*dir2);
	MGStraight* sl= new MGStraight();
	sl->set_straight(MGSTRAIGHT_UNLIMIT, *dir, origin);
	return sl;
}

//Compute part of the surface limitted by the parameter range bx.
//bx(0) is the parameter (us,ue) and bx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
MGSurface* MGPlane::part(const MGBox& uvbx, int multiple) const{
		MGSPointSeq sp(2, 2, 3);
		double u0=uvbx[0]. low_point(), u1=uvbx[0].high_point();
		double v0=uvbx[1]. low_point(), v1=uvbx[1].high_point();
		sp.store_at(0, 0, eval(u0,v0).data()); // u0 v0
		sp.store_at(1, 0, eval(u1,v0).data()); // u1 v0
		sp.store_at(0, 1, eval(u0,v1).data()); // u0 v1
		sp.store_at(1, 1, eval(u1,v1).data()); // u1 v1

		MGKnotVector tu(2, 2), tv(2, 2);
		tu[0] = tu[1] = u0;
		tu[2] = tu[3] = u1;
		tv[0] = tv[1] = v0;
		tv[2] = tv[3] = v1;
		return new MGSBRep(sp, tu, tv);
}

// Construct perimeter (u,v) parameter position.
// i is perimeter number:
// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
// t is perimeter parameter line's parameter value of u or v.
MGPosition MGPlane::perimeter_uv(unsigned i,double t) const{
	double u,v;
	double infinite;
	if(i==0 || i==3) infinite=-mgInfiniteVal;
	else infinite=mgInfiniteVal;
	if(i%2){
		u=t; v=infinite;
	}else{
		u=infinite; v=t;
	}
	return MGPosition(u,v);
}

// 与えられた点にもっとも近い面上の点とパラメータ値を返却する。
// Function's return value is if point is obtained(1) or not(0)
int MGPlane::perp_point (
	const MGPosition & p,			// 与えられた点
	MGPosition& uv,					// パラメータ値
	const MGPosition* uvguess 		// guess parameter value of the surface.
	)const{
	uv=param(p);
	return 1;
}

//Return all(actually one) foots of perpendicular straight lines from P.
MGPosition_list MGPlane::perps(
	const MGPosition& P				// Point of a space(指定点)
)const {
	MGPosition_list list;
	list.append(*this,param(P));
	return list;
}

//Test if the surface is planar or not.
//Returned is 0(false) if this is not planar, 1(true) if this planar.
int MGPlane::planar(
	MGPlane& plane,		//Plane that might be closest to this.
						//Plane is always output even if not planar.
	double& deviation	//maximum deviation of this from the output plane.
	) const{
	plane=*this;
	deviation=0.;
	return 1;
}

//Test if part of the surface is planar or not within the tolerance tol.
//The part of the surface is input by the surface parameter range uvbox.
//Returned is 0(false) if this is not planar, 1(true) if planar.
//For plane, planar always returns true.
int MGPlane::planar(
	const MGBox& uvbox,//This surface parameter range.
	double tol,	//maximum deviation allowed to regard the sub surface as a plane.
	int* divideU//Direction to subdivide will be output, if this was not planar.
				//=1: u direction, =0: v direction.
				) const{
	double u0=uvbox[0][0].value(), u1=uvbox[0][1].value();
	double v0=uvbox[1][0].value(), v1=uvbox[1][1].value();
	double um=(u0+u1)*0.5, vm=(v0+v1)*0.5;
	MGVector dfdu=eval(um,vm,1,0), dfdv=eval(um,vm,0,1);
	double ulen=dfdu.len()*(u1-u0), vlen=dfdv.len()*(v1-v0);

	if(divideU){
		if(ulen<vlen) *divideU=0;
		else *divideU=1;
	}
	return 1;

}

// 自身と与えられた平面との関係を返す。
MGPSRELATION MGPlane::relation(
	const MGPlane& pl,
	MGStraight& sl
	) const{
	if(m_normal==pl.m_normal){
		sl=MGStraight();
		if(MGAZero(m_d-pl.m_d)) return MGPSREL_COIN;
		else                    return MGPSREL_PARALLEL;
	}
	double g1[4], g2[4], line[2][3]; int iflag;
	for(size_t i=0; i<3; i++){
		g1[i]=m_normal.ref(i); g2[i]=pl.m_normal.ref(i);
	}
	g1[3]=m_d; g2[3]=pl.m_d;
	bvi2pl_(g1,g2,line[0],&iflag);	//Compute intersection line.
	MGVector dir(3,line[1]); MGPosition p(3,line[0]);
	sl=MGStraight(MGSTRAIGHT_UNLIMIT, dir, p);
	return MGPSREL_ISECT;
}

// 自身と与えられた直線の関係を返す。
MGPSRELATION MGPlane::relation (
	const MGStraight& s,
	MGCSisect& inter
	) const {
	return s.relation(*this,inter);
}

//Return the space dimension.
size_t MGPlane::sdim() const{
	size_t dim1=m_root_point.sdim(), dim2=m_uderiv.sdim();
	if(dim1<dim2) dim1=dim2;
	dim2=m_vderiv.sdim();
	if(dim1<dim2) dim1=dim2;
	return dim1;
}

//Obtain boundary and main parameter lines of the FSurface.
//skeleton includes boundary() and inner parameter lines.
//density indicates how many inner parameter lines are necessary
//for both u and v directions.
MGPvector<MGCurve> MGPlane::skeleton(int density)const{
	MGPvector<MGCurve> crv_list;
	crv_list.push_back(parameter_curve(true, 0.0));
	crv_list.push_back(parameter_curve(false, 0.0));
	return crv_list;
}

//Obtain all the parameter curves at knots of u and v knot vector.
MGPvector<MGCurve> MGPlane::skeleton_at_knots()const{
	MGPvector<MGCurve> crv_list;
	crv_list.push_back(parameter_curve(true, 0.0));
	crv_list.push_back(parameter_curve(false, 0.0));
	return crv_list;
}

//split this fsurface at the parameter param.
void MGPlane::split(
	double param,//parameter value of this fsurface. if is_u is true, param is u-value,
				//else v-value.
	bool is_u,	//indicates if param is u or v of the surface parameter (u,v).
	MGPvector<MGFSurface>& surfaces//splitted surfaces will be output.
)const{
	surfaces.clear();
	surfaces.push_back(clone());
}

// 点を平面に投影した点の平面のパラメータ表現(u,v)を求める。
MGPosition  MGPlane::uv(const MGPosition& p) const{
	MGVector vp(p-m_root_point);
	return uv(vp);
}

// Vectorを平面に投影したVectorの平面のパラメータ表現(u,v)を求める。
MGVector  MGPlane::uv(const MGVector& vec) const{
	double u,v;
	if(m_uderiv.orthogonal(m_vderiv)){
		u=(vec%m_uderiv)/(m_uderiv%m_uderiv);
		v=(vec%m_vderiv)/(m_vderiv%m_vderiv);
	}else{
		MGUnit_vector uunit,vunit;
		m_normal.orthonormal(m_uderiv,uunit,vunit);
		v=(vec%vunit)/(m_vderiv%vunit);
		m_normal.orthonormal(m_vderiv,vunit,uunit);
		u=(vec%uunit)/(m_uderiv%uunit);
	}
	return MGVector(u,v);
}

//////////Operator overload 演算子の多重定義/////////////

//Assignment.
//When the leaf object of this and srf2 are not equal, this assignment
//does nothing.
MGPlane& MGPlane::operator=(const MGPlane& pl){
	if(this==&pl)
		return *this;

	MGSurface::operator=(pl);
	if(m_uknotV)
		delete m_uknotV; m_uknotV=0;
	if(m_vknotV)
		delete m_vknotV; m_vknotV=0;
	m_d=pl.m_d;
	m_normal=pl.m_normal;
	m_root_point=pl.m_root_point;
	m_uderiv=pl.m_uderiv;
	m_vderiv=pl.m_vderiv;

	return *this;
}
MGPlane& MGPlane::operator=(const MGGel& gel2){
	const MGPlane* gel2_is_this=dynamic_cast<const MGPlane*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

// 平面に平行移動を行ないオブジェクトを生成する。
MGPlane MGPlane::operator+ (const MGVector& v) const{
	MGPlane pl = *this;
	pl += v;
	return pl;
}
MGPlane operator+ (const MGVector& v, const MGPlane& pl){
	return pl+v;
}

// 平面に平行移動を行ない自身の平面とする。
MGPlane& MGPlane::operator+= (const MGVector& v){
	m_root_point += v;
	m_d = m_normal%m_root_point;
	if(m_box) (*m_box)+=v;
	return *this;
}

// 平面に逆方向の平行移動を行ないオブジェクトを生成する。
MGPlane MGPlane::operator- (const MGVector& v) const{
	MGPlane pl = *this;
	pl -= v;
	return pl;
}

// 平面に逆方向の平行移動を行ない自身の平面とする。
MGPlane& MGPlane::operator-= (const MGVector& v) {
	m_root_point -= v;
	m_d = m_normal%m_root_point;
	if(m_box) (*m_box)-=v;
	return *this;
}

//平面のスケーリングを行い，平面を作成する。
//Scaling of the plane by a double.
MGPlane MGPlane::operator* (double scale) const{
	MGPlane pl=*this;
	pl *= scale;
	return pl;
}

//平面のスケーリングを行い，平面を作成する。
//Scaling of the plane by a double.
MGPlane operator* (double scale, const MGPlane& pl){
	return pl*scale;
}

//平面のスケーリングを行い自身の平面とする。
//Scaling of the plane by a double.
MGPlane& MGPlane::operator*= (double scale){
	m_root_point *= scale;
	m_uderiv*=scale; m_vderiv*=scale;
	m_d *= scale;
	update_mark();
	return *this;
}

// 与えられた変換で平面の変換を行い,平面を作成する。
MGPlane MGPlane::operator* (const MGMatrix& mat) const{
	MGPlane pl=*this;
	pl *= mat;
	return pl;
}

// 与えられた変換で平面の変換を行い自身の平面とする。
MGPlane& MGPlane::operator*= (const MGMatrix& mat){
	m_root_point *= mat;
	m_uderiv = m_uderiv * mat;
	m_vderiv = m_vderiv * mat;
	MGUnit_vector uunit=m_uderiv,vunit;
	uunit.orthonormal(m_vderiv, vunit, m_normal);
	m_d=m_root_point%m_normal;
	update_mark();
	return *this;
}

// 与えられた変換によってトランスフォームをおこない平面を生成する。
MGPlane MGPlane::operator* (const MGTransf& t) const{
	MGPlane pl = *this;
	pl *= t;
	return pl;
}

// 与えられた変換によってトランスフォームをおこない自身の平面にする。
MGPlane& MGPlane::operator*= (const MGTransf& t) {
	m_root_point *= t;
	m_uderiv = m_uderiv * t.affine();
	m_vderiv = m_vderiv * t.affine();
	MGUnit_vector uunit=m_uderiv, vunit;
	uunit.orthonormal(m_vderiv, vunit, m_normal);
	m_d=m_root_point%m_normal;
	update_mark();
	return *this;
}

// 論理演算子の多重定義
// 自身の平面と与えられた平面が等しいかどうか比較し判定する。
bool MGPlane::operator==(const MGPlane& srf2) const{
	if(m_normal.sdim()==0 && srf2.m_normal.sdim()==0)
		return true;
	return (m_normal==srf2.m_normal && MGAZero(m_d-srf2.m_d)) ;
}

bool MGPlane::operator<(const MGPlane& gel2)const{
	return m_d<gel2.m_d;
}
bool MGPlane::operator==(const MGGel& gel2)const{
	const MGPlane* gel2_is_this=dynamic_cast<const MGPlane*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGPlane::operator<(const MGGel& gel2)const{
	const MGPlane* gel2_is_this=dynamic_cast<const MGPlane*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}
