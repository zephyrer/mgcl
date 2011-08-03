/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGBSumSurf.cpp
//
// Defines Boolean sum surface.
//Boolian sum surface is defined by three surfaces g1, g2, and g12 as
//f(u,v)=g1+g2-g12. Typically Gordon surface is a boolean sum surface.
//See "Curves and Srufaces fro CAGD" by Gerald Farin.

////////// Constructor //////////

//construct from newed MGSurface. The ownership of g1, g2, and g12 will be
//transfered to this MGBSumSurf.
MGBSumSurf::MGBSumSurf(
	MGSurface* g1,
	MGSurface* g2,
	MGSurface* g12
):MGSurface(*g1),m_g1(g1),m_g2(g2),m_g12(g12){
	update_mark();
	assert(g1->param_range()==g2->param_range());
	assert(g2->param_range()==g12->param_range());
}

//construct from three MGSurface.
MGBSumSurf::MGBSumSurf(
	const MGSurface& g1,
	const MGSurface& g2,
	const MGSurface& g12
):MGSurface(g1){
	update_mark();
	assert(g1.param_range()==g2.param_range());
	assert(g2.param_range()==g12.param_range());
	m_g1=g1.copy_surface();
	m_g2=g2.copy_surface();
	m_g12=g12.copy_surface();
}

//Copy constructor.
MGBSumSurf::MGBSumSurf(const MGBSumSurf& rhs):MGSurface(rhs)
,m_g1(0),m_g2(0),m_g12(0){
	if(!rhs.m_g1) return;
	m_g1=rhs.m_g1->copy_surface();
	m_g2=rhs.m_g2->copy_surface();
	m_g12=rhs.m_g12->copy_surface();
}

//////////Destructor////////

MGBSumSurf::~MGBSumSurf(){
	delete m_g1;
	delete m_g2;
	delete m_g12;
}
									
////////// Operator Overload //////////

MGBSumSurf& MGBSumSurf::operator=(const MGBSumSurf& rhs){
	if(&rhs==this)
		return *this;

	MGSurface::operator=(rhs);
	delete m_g1; m_g1=0;
	delete m_g2; m_g2=0;
	delete m_g12; m_g12=0;

	if(!rhs.m_g1)
		return *this;
	m_g1=rhs.m_g1->copy_surface();
	m_g2=rhs.m_g2->copy_surface();
	m_g12=rhs.m_g12->copy_surface();
	return *this;
}
MGBSumSurf& MGBSumSurf::operator=(const MGGel& gel2){
	const MGBSumSurf* gel2_is_this=dynamic_cast<const MGBSumSurf*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

// 曲面の平行移動を行ない自身とする。
//Translation of the surface.
MGBSumSurf& MGBSumSurf::operator+= (const MGVector& v){
	(*m_g1)+=v;
	(*m_g2)+=v;
	(*m_g12)+=v;
	return *this;
}

// 面の逆方向の平行移動を行ない自身の面とする。
//Translation of the surface.
MGBSumSurf& MGBSumSurf::operator-= (const MGVector& v){
	(*m_g1)-=v;
	(*m_g2)-=v;
	(*m_g12)-=v;
	return *this;
}
												  
// 与えられたスケーリングで曲面の変換を行い自身の曲面とする。
//Scaling of the surface.
MGBSumSurf& MGBSumSurf::operator*= (double scale){
	(*m_g1)*=scale;
	(*m_g2)*=scale;
	(*m_g12)*=scale;
	return *this;
}
												  
// 与えられた変換で曲面の変換を行い自身の曲面とする。
//Matrix transformation of the surface.
MGBSumSurf& MGBSumSurf::operator*= (const MGMatrix& mat){
	(*m_g1)*=mat;
	(*m_g2)*=mat;
	(*m_g12)*=mat;
	return *this;
}

// 与えられた変換による面のトランスフォームを行ない自身とする。
//Matrix transformation of the surface.
MGBSumSurf& MGBSumSurf::operator*= (const MGTransf& tr){
	(*m_g1)*=tr;
	(*m_g2)*=tr;
	(*m_g12)*=tr;
	return *this;
}

//Equal operator overload. 論理演算子多重定義
//Comparison of two surfaces. 自身と与曲面が等しいか比較する。
bool MGBSumSurf::operator==(const MGBSumSurf& srf2)const{
	if(!m_g1 && !srf2.m_g1) return true;
	if((*m_g1)!=(*(srf2.m_g1))) return false;
	if((*m_g2)!=(*(srf2.m_g2))) return false;
	if((*m_g12)!=(*(srf2.m_g12))) return false;
	return true;
}
bool MGBSumSurf::operator<(const MGBSumSurf& gel2)const{
	return (*m_g1)<*(gel2.m_g1);
}
bool MGBSumSurf::operator==(const MGGel& gel2)const{
	const MGBSumSurf* gel2_is_this=dynamic_cast<const MGBSumSurf*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGBSumSurf::operator<(const MGGel& gel2)const{
	const MGBSumSurf* gel2_is_this=dynamic_cast<const MGBSumSurf*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

////////// Member Function //////////

size_t MGBSumSurf::bdim_u() const	//Returns B-Rep Dimension of u.
{
	const MGKnotVector& t=knot_vector_u();
	return t.bdim();
}
size_t MGBSumSurf::bdim_v() const	//Returns B-Rep Dimension of v.
{
	const MGKnotVector& t=knot_vector_v();
	return t.bdim();
}

// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
//Return minimum box that includes limitted surface by uvrange.
MGBox MGBSumSurf::box_limitted(
	const MGBox& uvrange	// Parameter Range of the curve.
)const{
	MGBox bx=m_g1->box_limitted(uvrange);
	bx|=m_g2->box_limitted(uvrange);
	bx|=m_g12->box_limitted(uvrange);
	return bx;
}

//Changing this object's space dimension.
MGBSumSurf& MGBSumSurf::change_dimension(
	size_t sdim,	// new space dimension
	size_t start1, 	// Destination order of new object.
	size_t start2)	// Source order of this object.
{
	m_g1->change_dimension(sdim,start1,start2);
	m_g2->change_dimension(sdim,start1,start2);
	m_g12->change_dimension(sdim,start1,start2);
	return *this;
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
MGBSumSurf& MGBSumSurf::change_range(		//BLUCPR
	int is_u,				//if true, (t1,t2) are u-value. if not, v.
	double t1,				//Parameter value for the start of original. 
	double t2				//Parameter value for the end of original. 
){
	m_g1->change_range(is_u,t1,t2);
	m_g2->change_range(is_u,t1,t2);
	m_g12->change_range(is_u,t1,t2);
	return *this;
}

//Return minimum box that includes whole of the surface.
//Returned is a newed object pointer.
MGBox* MGBSumSurf::compute_box()const{
	MGBox* bx=new MGBox;
	(*bx)|=m_g1->box();
	(*bx)|=m_g2->box();
	(*bx)|=m_g12->box();
	return bx;
}

//Construct new surface object by copying to newed area.
//User must delete this copied object by "delete".
MGBSumSurf* MGBSumSurf::clone() const{
	MGBSumSurf* bss=new MGBSumSurf(*this);
	return bss;
}

//Construct new surface object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGBSumSurf* MGBSumSurf::copy_change_dimension(
	size_t sdim,			// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2 		// Source order of this line.
)const{
	MGBSumSurf* bss=new MGBSumSurf(*this);
	bss->change_dimension(sdim,start1,start2);
	return bss;
}

//Evaluate surface data.
//Currently ndu=ndv=0 is assumed.
MGVector MGBSumSurf::eval(
	double u, double v	// Parameter value of the surface.
						// must be 0<=u,v<=1.
	, size_t ndu		// Order of derivative along u.
	, size_t ndv		// Order of derivative along v.
)const{
	MGVector val1=m_g1->eval(u,v,ndu,ndv);
	MGVector val2=m_g2->eval(u,v,ndu,ndv);
	MGVector val12=m_g12->eval(u,v,ndu,ndv);
	//cout<<val1<<val2<<val12<<endl;
	return val1+val2-val12;
}

// Exchange parameter u and v.
MGSurface& MGBSumSurf::exchange_uv(){
	m_g1->exchange_uv();
	m_g2->exchange_uv();
	m_g12->exchange_uv();
	return *this;
}

bool MGBSumSurf::in_range(double u, double v) const{
	return m_g1->in_range(u,v);
}

//The following two function will be used in perps or isect
//to decide how many division of the surface along u or v direction
//should be applied before using perp_guess or isect_guess.
size_t MGBSumSurf::intersect_dnum_u() const{
	size_t n=m_g1->intersect_dnum_u();
	size_t n2=m_g2->intersect_dnum_u();
	if(n<n2) n=n2;
	size_t n3=m_g12->intersect_dnum_u();
	if(n<n3) n=n3;
	return n;
}
size_t MGBSumSurf::intersect_dnum_v() const{
	size_t n=m_g1->intersect_dnum_v();
	size_t n2=m_g2->intersect_dnum_v();
	if(n<n2) n=n2;
	size_t n3=m_g12->intersect_dnum_v();
	if(n<n3) n=n3;
	return n;
}

MGCSisect_list MGBSumSurf::isect(const MGCurve& curve)const{
	return curve.isect(*this);
}

MGCSisect_list MGBSumSurf::isect(const MGRLBRep& curve)const{
	return curve.isect(*this);
}

MGCSisect_list MGBSumSurf::isect(const MGEllipse& curve)const{
	return curve.isect(*this);
}

MGCSisect_list MGBSumSurf::isect(const MGLBRep& curve)const{
	return curve.isect(*this);
}

MGCSisect_list MGBSumSurf::isect(const MGSurfCurve& curve)const{
	return curve.isect(*this);
}

MGCSisect_list MGBSumSurf::isect(const MGBSumCurve& curve)const{
	return curve.isect(*this);
}

// Surface と Surface の交線を求める。
//Compute intesection of Sphere and Surface.
MGSSisect_list MGBSumSurf::isect(const MGSurface& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGBSumSurf::isect(const MGPlane& srf2) const{
	return intersectPl(srf2);
}
MGSSisect_list MGBSumSurf::isect(const MGSphere& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGBSumSurf::isect(const MGCylinder& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGBSumSurf::isect(const MGSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGBSumSurf::isect(const MGRSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGBSumSurf::isect(const MGBSumSurf& srf2) const{
	return intersect(srf2);
}

//"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGBSumSurf::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	size_t incr			//Incremental valuse of B-coef's id.
)const{
	MGCurve* g1=m_g1->isect_incr_pline(uv,kdt,du,dv,u,v,incr);
	MGCurve* g2=m_g2->isect_incr_pline(uv,kdt,du,dv,u,v,incr);
	MGCurve* g12=m_g12->isect_incr_pline(uv,kdt,du,dv,u,v,incr);
	return new MGBSumCurve(g1,g2,g12);
}

//Return order of intersection line order of MGLBRep.
//The default is 4.
size_t MGBSumSurf::isect_order() const{
	size_t n=m_g1->isect_order();
	size_t n2=m_g2->isect_order();
	if(n<n2) n=n2;
	size_t n3=m_g12->isect_order();
	if(n<n3) n=n3;
	return n;
}

//Access to i-th element of u knot
double MGBSumSurf::knot_u(size_t i) const{
	const MGKnotVector& t=knot_vector_u();
	return t[i];
}

//Access to i-th element of v knot
double MGBSumSurf::knot_v(size_t i) const{
	const MGKnotVector& t=knot_vector_v();
	return t[i];
}

//Returns the u knot vector.
const MGKnotVector& MGBSumSurf::knot_vector_u()const{
	size_t nu1=(m_g1->knot_vector_u()).bdim();
	size_t nu2=(m_g2->knot_vector_u()).bdim();
	size_t nu3=(m_g12->knot_vector_u()).bdim();
	if(nu1>=nu2){
		if(nu1>=nu3) return m_g1->knot_vector_u();
		else return m_g12->knot_vector_u();
	}else if(nu2>=nu3) return m_g2->knot_vector_u();
	else return m_g12->knot_vector_u();
}

MGKnotVector& MGBSumSurf::knot_vector_u(){
	const MGBSumSurf*bss=const_cast<const MGBSumSurf*>(this);
	const MGKnotVector& t=bss->knot_vector_u();
	return *(const_cast<MGKnotVector*>(&t));
}

//Returns the v knot vector.
const MGKnotVector& MGBSumSurf::knot_vector_v() const{
	size_t nu1=(m_g1->knot_vector_v()).bdim();
	size_t nu2=(m_g2->knot_vector_v()).bdim();
	size_t nu3=(m_g12->knot_vector_v()).bdim();
	if(nu1>=nu2){
		if(nu1>=nu3) return m_g1->knot_vector_v();
		else return m_g12->knot_vector_v();
	}else if(nu2>=nu3) return m_g2->knot_vector_v();
	else return m_g12->knot_vector_v();
}
MGKnotVector& MGBSumSurf::knot_vector_v(){
	const MGBSumSurf*bss=const_cast<const MGBSumSurf*>(this);
	const MGKnotVector& t=bss->knot_vector_v();
	return *(const_cast<MGKnotVector*>(&t));
}

//Negate direction of surface.
void MGBSumSurf::negate(
	int is_u)	// Negate along u-direction if is_u is ture,
				// else along v-direction.
{
	m_g1->negate(is_u);
	m_g2->negate(is_u);
	m_g12->negate(is_u);
}

unsigned MGBSumSurf::order_u() const	//Returns the order of u.
{
	const MGKnotVector& t=knot_vector_u();
	return t.order();
}
unsigned MGBSumSurf::order_v() const	//Returns the order of v.
{
	const MGKnotVector& t=knot_vector_v();
	return t.order();
}

// Compute parameter curve.
//Returned is newed area pointer, and must be freed by delete.
MGCurve* MGBSumSurf::parameter_curve(
	int is_u				//Indicates x is u-value if is_u is true.
	, double x				//Parameter value.
							//The value is u or v according to is_u.
)const{
	MGCurve* g1=m_g1->parameter_curve(is_u,x);
	MGCurve* g2=m_g2->parameter_curve(is_u,x);
	MGCurve* g12=m_g12->parameter_curve(is_u,x);
	return new MGBSumCurve(g1,g2,g12);
}

// Return ending parameter value.
double MGBSumSurf::param_e_u()const{
	return m_g1->param_e_u();
}

double MGBSumSurf::param_e_v() const{
	return m_g1->param_e_v();
}

// パラメータ範囲を返す。
//Return parameter range.
MGBox MGBSumSurf::param_range() const{
	return m_g1->param_range();
}

// Return starting parameter value.
double MGBSumSurf::param_s_u() const{
	return m_g1->param_s_u();
}
double MGBSumSurf::param_s_v() const{
	return m_g1->param_s_v();
}

//Compute part of the surface limitted by the parameter range bx.
//bx(0) is the parameter (us,vs) and bx(1) is (ue,ve).
//That is u range is from us to ue , and so on.
//Retured is newed object, must be deleted.
MGBSumSurf* MGBSumSurf::part(
	const MGBox& bx,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
)const{
	MGSurface* g1=m_g1->part(bx,multiple);
	MGSurface* g2=m_g2->part(bx,multiple);
	MGSurface* g12=m_g12->part(bx,multiple);
	return new MGBSumSurf(g1,g2,g12);
}

//Retrieve perimeter i of this surface.
// i must be < perimeter_num().
//When perimeter_num()==0, this function is undefined.
//Retured is newed object, must be deleted.
MGCurve* MGBSumSurf::perimeter_curve(size_t i)const{
	assert(i<4);

	int is_u; double x;
	switch(i){
		case 0:  is_u=0; x=param_s_v(); break;
		case 1:  is_u=1; x=param_e_u(); break;
		case 2:  is_u=0; x=param_e_v(); break;
		default: is_u=1; x=param_s_u(); break;
	}
	return parameter_curve(is_u, x);
}

//Shrink this surface to the part limitted by the parameter range of uvbx.
//New parameter range uvbx2 is so determined that uvbx2 is the smallest
//box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
//the values of u or v knots of the surface knotvector.
//uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
void MGBSumSurf::shrink_to_knot(
	const MGBox& uvbx,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
){
	m_g1->shrink_to_knot(uvbx,multiple);
	m_g2->shrink_to_knot(uvbx,multiple);
	m_g12->shrink_to_knot(uvbx,multiple);
}

//Obtain 1D surface rep. of this surf which can be used for
//isect(const MGPlane& pl). This surf1D is used in isect for
//the argument of isect_startPlane, which will use surf1D to compute isect(pl).
//surf1D=0.(intersection with x=0. plane) is the intersection lines.
MGSBRep* MGBSumSurf::surf1D(
	const MGPlane& pl
)const{
	assert(false); return 0;
}

//メンバデータを読み込む関数
void MGBSumSurf::ReadMembers(MGIfstream& buf){
	MGSurface::ReadMembers(buf);
	m_g1=static_cast<MGBSumSurf*>(buf.ReadPointer());
	m_g2=static_cast<MGBSumSurf*>(buf.ReadPointer());
	m_g12=static_cast<MGBSumSurf*>(buf.ReadPointer());
}
	
//メンバデータを書き込む関数
void MGBSumSurf::WriteMembers(MGOfstream& buf) const{
	MGSurface::WriteMembers(buf);
	buf.WritePointer(m_g1);
	buf.WritePointer(m_g2);
	buf.WritePointer(m_g12);
}
