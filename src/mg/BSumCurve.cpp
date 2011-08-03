/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#include "mg/BSumCurve.h"
#include "mg/LBRep.h"
#include "mg/Straight.h"
#include "mg/RLBRep.h"
#include "mg/Ellipse.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/Surface.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Define MGBSumCurve Class(Boolean sum curve of three curves).
//MGBSumCurve is a curve to support the class MGBSumSurf, is boolean sum of
//three curves, m_g1, m_g2, and m_g12.
//MGBSumCurve(t) is defined as m_g1(t)+m_g2(t)-m_g12(t).
//

////////////Constructor/////////////

//Void constructor(初期化なしでオブジェクトを作成する。)
MGBSumCurve::MGBSumCurve():MGCurve(),m_g1(0),m_g2(0),m_g12(0){;}

//Copy constructor.
MGBSumCurve::MGBSumCurve(const MGBSumCurve& curve)
:MGCurve(curve),m_g1(0),m_g2(0),m_g12(0){
	if(!curve.m_g1) return;
	m_g1=curve.m_g1->clone();
	m_g2=curve.m_g2->clone();
	m_g12=curve.m_g12->clone();
}

//constructor of three curves.
//The ownership of g1, g2, and g12 will be transfered to this MGBSumCurve.
MGBSumCurve::MGBSumCurve(MGCurve* g1, MGCurve* g2, MGCurve* g12)
:MGCurve(*g1),m_g1(g1), m_g2(g2), m_g12(g12){
	update_mark();
}

//constructor of three curves.
//The ownership of g1, g2, and g12 will be transfered to this MGBSumCurve.
MGBSumCurve::MGBSumCurve(const MGCurve& g1, const MGCurve& g2, const MGCurve& g12)
:MGCurve(g1),m_g1(g1.clone()), m_g2(g2.clone()), m_g12(g12.clone()){
	update_mark();
}

////////// Destructor //////////
MGBSumCurve::~MGBSumCurve(){
	delete m_g1;
	delete m_g2;
	delete m_g12;
}

////////// Operator overload(演算子多重定義) //////////

//Assignment.
//When the leaf object of this and geo2 are not equal, this assignment
//does nothing.
MGBSumCurve& MGBSumCurve::operator=(const MGBSumCurve& crv2){
	if(&crv2==this)
		return *this;

	MGCurve::operator=(crv2);
	delete m_g1; m_g1=0;
	delete m_g2; m_g2=0;
	delete m_g12; m_g12=0;
	if(!crv2.m_g1) return *this;
	m_g1=crv2.m_g1->clone();
	m_g2=crv2.m_g2->clone();
	m_g12=crv2.m_g12->clone();
	return *this;
}
MGBSumCurve& MGBSumCurve::operator=(const MGGel& gel2){
	const MGBSumCurve* gel2_is_this=dynamic_cast<const MGBSumCurve*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線を平行移動して自身とする。
MGBSumCurve& MGBSumCurve::operator+= (const MGVector& v){
	(*m_g1)+=v;
	(*m_g2)+=v;
	(*m_g12)+=v;
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGBSumCurve& MGBSumCurve::operator-= (const MGVector& v){
	(*m_g1)-=v;
	(*m_g2)-=v;
	(*m_g12)-=v;
	return *this;
}

//Update the curve by multiplying scale.
// 与えられたスケールを曲線にかける。
MGBSumCurve& MGBSumCurve::operator*= (double scale){
	(*m_g1)*=scale;
	(*m_g2)*=scale;
	(*m_g12)*=scale;
	return *this;
}

//Update the curve by transformation of matrix.
// 与えられた変換で直線の変換を行い自身の直線とする。
MGBSumCurve& MGBSumCurve::operator*= (const MGMatrix& mat){
	(*m_g1)*=mat;
	(*m_g2)*=mat;
	(*m_g12)*=mat;
	return *this;
}

//Update the curve by transformation of transf.
// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGBSumCurve& MGBSumCurve::operator*= (const MGTransf& tr){
	(*m_g1)*=tr;
	(*m_g2)*=tr;
	(*m_g12)*=tr;
	return *this;
}

//Transformation object construction
MGBSumCurve MGBSumCurve::operator+ (const MGVector& v)const{
	MGBSumCurve cc(*this);
	cc+=v;
	return cc;
}
MGBSumCurve operator+ (const MGVector& v, const MGBSumCurve& cc){
	return cc+v;
}
MGBSumCurve MGBSumCurve::operator- (const MGVector& v) const{
	MGBSumCurve cc(*this);
	cc-=v;
	return cc;
}
MGBSumCurve MGBSumCurve::operator* (double scale) const{
	MGBSumCurve cc(*this);
	cc*=scale;
	return cc;
}
MGBSumCurve operator* (double scale, const MGBSumCurve& cc){
	return cc*scale;
}
MGBSumCurve MGBSumCurve::operator* (const MGMatrix& mat) const{
	MGBSumCurve cc(*this);
	cc*=mat;
	return cc;
}
MGBSumCurve MGBSumCurve::operator* (const MGTransf& tr) const{
	MGBSumCurve cc(*this);
	cc*=tr;
	return cc;
}

//Logical operator overload(論理演算子多重定義)
//Test if two curves are equal.
// 与曲線と自身が等しいかの比較判定を行う。
bool MGBSumCurve::operator== (const MGBSumCurve& crv2)const{
	if(!m_g1 && !crv2.m_g1)
		return true;
	if((*m_g1)!=(*(crv2.m_g1)))
		return false;
	if((*m_g2)!=(*(crv2.m_g2)))
		return false;
	if((*m_g12)!=(*(crv2.m_g12)))
		return false;
	return true;
}
bool MGBSumCurve::operator<(const MGBSumCurve& gel2)const{
	return *m_g1<*(gel2.m_g1);
}
bool MGBSumCurve::operator==(const MGGel& gel2)const{
	const MGBSumCurve* gel2_is_this=dynamic_cast<const MGBSumCurve*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGBSumCurve::operator<(const MGGel& gel2)const{
	const MGBSumCurve* gel2_is_this=dynamic_cast<const MGBSumCurve*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

////////// Member Function ///////////

//Returns B-Rep Dimension.
size_t MGBSumCurve::bdim()const{
	const MGKnotVector& t=knot_vector();
	return t.bdim();
}

//Return minimum box that includes the curve of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox MGBSumCurve::box_limitted(
	const MGInterval& rang // Parameter Range of the curve.
)const{
	MGBox bx=m_g1->box_limitted(rang);
	bx|=m_g2->box_limitted(rang);
	bx|=m_g12->box_limitted(rang);
	return bx;
}

//Changing this object's space dimension.
MGBSumCurve& MGBSumCurve::change_dimension(
	size_t sdim,	// new space dimension
	size_t start1, 	// Destination order of new object.
	size_t start2 	// Source order of this object.
){
	m_g1->change_dimension(sdim,start1,start2);
	m_g2->change_dimension(sdim,start1,start2);
	m_g12->change_dimension(sdim,start1,start2);
	return *this;
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
void MGBSumCurve::change_range(
	double t1,	//Parameter value for the start of original. 
	double t2	//Parameter value for the end of original. 
){
	m_g1->change_range(t1,t2);
	m_g2->change_range(t1,t2);
	m_g12->change_range(t1,t2);
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGBSumCurve& MGBSumCurve::coordinate_exchange(
	size_t i, size_t j
){
	m_g1->coordinate_exchange(i,j);
	m_g2->coordinate_exchange(i,j);
	m_g12->coordinate_exchange(i,j);
	return *this;
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGBSumCurve* MGBSumCurve::clone()const{
	MGBSumCurve* bss=new MGBSumCurve(*this);
	return bss;
}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGBSumCurve::copy_as_nurbs()const{
	return new MGLBRep(*this);
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGBSumCurve* MGBSumCurve::copy_change_dimension(
	size_t sdim,	// new space dimension
	size_t start1,	// Destination order of new line.
	size_t start2 	// Source order of this line.
)const{
	MGBSumCurve* bss=new MGBSumCurve(*this);
	bss->change_dimension(sdim,start1,start2);
	return bss;
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGBSumCurve::eval(
	double t,		//Parameter value.
	size_t nderiv,	//Order of Derivative.
	int left		//Left continuous(left=true)
					//or right continuous(left=false).
)const{
	MGVector val=m_g1->eval(t,nderiv,left);
	val+=m_g2->eval(t,nderiv,left);
	val-=m_g12->eval(t,nderiv,left);
	return val;
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGBSumCurve::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){

}

//Intersection of MGBSumCurve and curve.
MGCCisect_list MGBSumCurve::isect(const MGCurve& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection of MGBSumCurve and MGStraight.
MGCCisect_list MGBSumCurve::isect(const MGStraight& curve2)const{
	return intersect(curve2);
}

//Intersection of MGBSumCurve and MGSurfCurve.
MGCCisect_list MGBSumCurve::isect(const MGSurfCurve& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

MGCSisect_list MGBSumCurve::isect(const MGSurface& srf)const{
	return srf.isect(*this);
}

//Intersection with a Surface
MGCSisect_list MGBSumCurve::isect(const MGPlane& surf)const{
	return intersect_with_plane(surf);
}

//Intersection with a Surface
MGCSisect_list MGBSumCurve::isect(const MGSphere& surf)const{
	return surf.intersect(*this);
}

//Intersection with a Surface
MGCSisect_list MGBSumCurve::isect(const MGCylinder& surf)const{
	return surf.intersect(*this);
}

//Intersection with a Surface
MGCSisect_list MGBSumCurve::isect(const MGSBRep& surf)const{
	return surf.intersect(*this);
}

//Intersection with a Surface
MGCSisect_list MGBSumCurve::isect(const MGRSBRep& surf)const{
	return surf.intersect(*this);
}

//Intersection with a Surface
MGCSisect_list MGBSumCurve::isect(const MGBSumSurf& surf)const{
	return surf.intersect(*this);
}

//Access to i-th element of knot
double MGBSumCurve::knot(size_t i)const{
	const MGKnotVector& t=knot_vector();
	return t[i];
}

//Returns the knot vector of the curve.
const MGKnotVector& MGBSumCurve::knot_vector() const{
	size_t nu1=(m_g1->knot_vector()).bdim();
	size_t nu2=(m_g2->knot_vector()).bdim();
	size_t nu3=(m_g12->knot_vector()).bdim();
	if(nu1>=nu2){
		if(nu1>=nu3) return m_g1->knot_vector();
		else return m_g12->knot_vector();
	}else if(nu2>=nu3) return m_g2->knot_vector();
	else return m_g12->knot_vector();
}

//Update this by limiting the parameter range of the curve.
// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
MGBSumCurve& MGBSumCurve::limit(const MGInterval& rng){
	m_g1->limit(rng);
	m_g2->limit(rng);
	m_g12->limit(rng);
	return *this;
}

//Negate the curve direction(曲線の方向を反転する)
void MGBSumCurve::negate(){
	m_g1->negate();
	m_g2->negate();
	m_g12->negate();
}

//Obtain parameter value if this curve is negated by "negate()".
double MGBSumCurve::negate_param(double t)const{
	double tspte=param_s()+param_e();
	return tspte-t;
}

//Returns the order.
unsigned MGBSumCurve::order()const{
	const MGKnotVector& kv=knot_vector();
	return kv.order();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGBSumCurve::param_normalize(double t)const{
	const MGKnotVector& kv=knot_vector();
	return kv.param_normalize(t);
}

//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and should be deleted
//by calling program, or memory leaked.
MGBSumCurve* MGBSumCurve::part(
	double t1, double t2,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
)const{
	MGCurve* g1=m_g1->part(t1,t2,multiple);
	MGCurve* g2=m_g2->part(t1,t2,multiple);
	MGCurve* g12=m_g12->part(t1,t2,multiple);
	return new MGBSumCurve(g1,g2,g12);
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list MGBSumCurve::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGBSumCurve::perps(
	const MGStraight& crv2		//The second curve
)const{
	return perpsSl(crv2);	
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
MGSurface* MGBSumCurve::sweep(
	const MGUnit_vector& uvec,	//Sweep Direction.
	double start_dist,			//distance to start edge.
	double end_dist			//distance to end edge.
)const{
	MGLBRep tempCrv(*this);
	MGSurface *rtnSrf = tempCrv.sweep(uvec, start_dist, end_dist);
	return rtnSrf;
}

//Unlimit parameter range of the curve(limitをはずす)
MGCurve& MGBSumCurve::unlimit(){
	m_g1->unlimit();
	m_g2->unlimit();
	m_g12->unlimit();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGBSumCurve::unlimit_end(){
	m_g1->unlimit_end();
	m_g2->unlimit_end();
	m_g12->unlimit_end();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGBSumCurve::unlimit_start(){
	m_g1->unlimit_start();
	m_g2->unlimit_start();
	m_g12->unlimit_start();
	return *this;
}

// Output function.
std::ostream& MGBSumCurve::out(std::ostream& ostrm) const{
	ostrm<<"MGBSumCurve::"<<this;
	MGCurve::out(ostrm);
	ostrm<<",m_g1="<<m_g1;if(m_g1) ostrm<<","<<(*m_g1);
	ostrm<<",m_g2="<<m_g2;if(m_g2) ostrm<<","<<(*m_g2);
	ostrm<<",m_g12="<<m_g12;if(m_g12) ostrm<<","<<(*m_g12);
	return ostrm;
}

//メンバデータを読み出す関数
// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void MGBSumCurve::ReadMembers(MGIfstream& buf){
	MGCurve::ReadMembers(buf);
	m_g1=static_cast<MGBSumCurve*>(buf.ReadPointer());
	m_g2=static_cast<MGBSumCurve*>(buf.ReadPointer());
	m_g12=static_cast<MGBSumCurve*>(buf.ReadPointer());
}

//メンバデータを書き込む関数
// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void MGBSumCurve::WriteMembers(MGOfstream& buf) const{
	MGCurve::WriteMembers(buf);
	buf.WritePointer(m_g1);
	buf.WritePointer(m_g2);
	buf.WritePointer(m_g12);
}

//Provide divide number of curve span for function intersect.
size_t MGBSumCurve::intersect_dnum() const{
	size_t n=m_g1->intersect_dnum();
	size_t n2=m_g2->intersect_dnum();
	if(n<n2) n=n2;
	size_t n3=m_g12->intersect_dnum();
	if(n<n3) n=n3;
	return n;
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGBSumCurve::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
)const{
	std::auto_ptr<MGCurve> g1=m_g1->oneD(g);
	std::auto_ptr<MGCurve> g2=m_g2->oneD(g);
	std::auto_ptr<MGCurve> g12=m_g12->oneD(g);
	return std::auto_ptr<MGCurve>(
		new MGBSumCurve(g1.release(),g2.release(),g12.release()));
}

//Return minimum box that includes whole of the curve.
//曲線部分を囲むボックスを返す。
MGBox* MGBSumCurve::compute_box() const{
	MGBox* bx=new MGBox;
	(*bx)|=m_g1->box();
	(*bx)|=m_g2->box();
	(*bx)|=m_g12->box();
	return bx;
}
