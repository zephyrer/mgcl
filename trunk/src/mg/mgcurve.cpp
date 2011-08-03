/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Unit_vector.h"
#include "mg/Transf.h"
#include "mg/Position.h"
#include "mg/Point.h"
#include "mg/Curve.h"
#include "mg/Ellipse.h"
#include "mg/PPRep.h"
#include "mg/RLBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CCisect_list.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/CompositeCurve.h"
#include "mg/Plane.h"
#include "mg/nlbit.h"
#include "mg/Tolerance.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGCurve
// Implementation of MGCurve.

// コンストラクタ
// 初期化なしでオブジェクトを作成する。
MGCurve::MGCurve():MGGeometry(){;}

//Copy constructor.
MGCurve::MGCurve(const MGCurve& curve):MGGeometry(curve){;}

// 仮想デストラクタ
MGCurve::~MGCurve(){;}

//Test if two curves are equal.
// 与曲線と自身が等しいかの比較判定を行う。
bool MGCurve::operator== (const MGCompositeCurve& crv) const{
	return crv.is_same_curve(*this);
}
bool MGCurve::operator== (const MGTrimmedCurve& crv) const{
	return crv.is_same_curve(*this);
}

//Approximate this curve as a MGLBRep curve from knot_vector[is] to [ie].
//This is an internal program of MGLBRep constructor.
void MGCurve::approximate_as_LBRep(
	MGLBRep& lb,		//Approximated obrep will be set.
	size_t order,		//new order
	size_t is, size_t ie//approximation parameter range, from lb.knot_vector()[is] to [ie].
)const{
	const MGKnotVector& tempKnotVector=knot_vector();
	size_t kold=tempKnotVector.order();
	size_t nnew=ie-is+kold-1;
	MGKnotVector ttmp(is-kold+1,nnew,tempKnotVector);
	MGKnotVector t(ttmp,order);//指定オーダーのノットベクトルに作り替える
	double err=t.param_error();
	double te=tempKnotVector(is);
	for(size_t i=is+1; i <=ie; i++){
		double	ts = te;			//スパンの始点
		te = tempKnotVector(i);		//スパンの終点
		if(te - ts <=err) continue;	//マルチノットのときの処理
		//1スパンの分割数を決定する
		int ndiv = offset_div_num(MGInterval(ts, te));//オフセット曲線の分割数
		double delta = (te - ts) / ndiv;
		for(int j=1; j<ndiv; j++){
			ts += delta;t.add_data(ts);
		}
	}

	//制御点を生成する
	MGNDDArray tau; tau.update_from_knot(t);
	size_t len = tau.length();
	MGBPointSeq bp1(len, sdim());
	for(size_t i = 0; i < len; i++)
		bp1.store_at(i, eval(tau(i)));

	//精度十分の曲線を生成する
	int error, fNRat = 0;
	const MGRLBRep *prbrep = dynamic_cast<const MGRLBRep*>(this);
	if(prbrep){
		if(prbrep->non_rational())
			fNRat=1;
	}

	//LBRepとRLBRepで非有理のときは、元のノットベクトルを使用する
	if(dynamic_cast<const MGLBRep*>(this) || fNRat){
		lb=MGLBRep(tau, bp1, t, error);
	}else{	//上記以外のとき点列から作り直す
		MGKnotVector tnew(tau,order);
		lb=MGLBRep(tau, bp1, tnew, error);
	}
	lb.remove_knot();
}

//Return minimum box that includes the geometry of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
//MGBox MGCurve::box(
//	const MGBox& bx	// Parameter Range of the geometry.
//					//bx's space dimension is manifold dimension
//					//of the geometry.
//) const{return box(bx(0));}

//Compute the point on the curve which is the intersection
//with x=point[0], or y=point[1] and the nearest point from point.
//If no intersection found, the nearer point parameter value out of start or end
//will be returned.
double MGCurve::closest2D(const MGPosition& point, double& dist)const{
	double x=point[0], y=point[1];
	MGCParam_list list=isect_1D(x);
	list.append(isect_1D(y,1));

	MGPosition Q;
	double t1,t,dist1;

	t=param_s();Q=eval(t);
	dist=fabs(Q[0]-x)+fabs(Q[1]-y);
	t1=param_e();Q=eval(t1);
	dist1=fabs(Q[0]-x)+fabs(Q[1]-y);
	if(dist1<dist) {t=t1; dist=dist1;}
	if(!list.size()) return t;

	MGCParam_list::Citerator i=list.begin(), ie=list.end();
	for(; i!=ie; i++){
		t1=*i;Q=eval(t1);
		dist1=fabs(Q[0]-x)+fabs(Q[1]-y);
		if(dist1<dist) {t=t1; dist=dist1;}
	}
	return t;
}

//Compute the closest point parameter value of this curve from a point.
double MGCurve::closest(const MGPosition& point) const{
	MGCParam_list list=perps(point); list.append(param_e());
	double t1,t=param_s();		//t is initial parameter value;
	double dist1,dist=(point-eval(t)).len();
	
	size_t n=list.entries();
	for(size_t i=0; i<n; i++){
		t1=list.removeFirst(); dist1=(point-eval(t1)).len();
		if(dist1<dist) {t=t1; dist=dist1;}
	}
	return t;
}

//Compute the closest point on the curve from this point.
//Function's return value is the parameter value of the curve.
double MGPosition::closest(const MGCurve& curve) const{
	return curve.closest(*this);
}

//Obtain ceter coordinate of the geometry.
MGPosition MGCurve::center() const{
	return eval((param_s()+param_e()) * 0.5);
}

//Obtain ceter parameter value of the geometry.
MGPosition MGCurve::center_param() const{
	double t=(param_s()+param_e())*0.5;
	return MGPosition(1,&t);
}

//Compute the closest point parameter value pair of this curve and crv2.
//MGPosition P of the function return contains this and crv2's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition MGCurve::closest(const MGCurve& crv2) const{
	MGPosition_list list=perps(crv2);
	double s0=param_s(), s1=param_e(), t0=crv2.param_s(), t1=crv2.param_e();
	MGPosition P0(eval(s0)), P1(eval(s1)), Q0(crv2.eval(t0)), Q1(crv2.eval(t1));
	list.append(MGPosition(closest(Q0), t0));
	list.append(MGPosition(closest(Q1), t1));
	list.append(MGPosition(s0,crv2.closest(P0)));
	MGPosition uv1;
	MGPosition uv(MGPosition(s1,crv2.closest(P1)));//uv is initial parameter value;
	double dist1,dist=(P1-crv2.eval(uv[1])).len();
	
	size_t n=list.entries();
	for(size_t i=0; i<n; i++){
		uv1=list.removeFirst();
		dist1=(eval(uv1.ref(0))-crv2.eval(uv1.ref(1))).len();
		if(dist1<dist) {uv=uv1; dist=dist1;}
	}
	return uv;
}

//曲線がCn連続かどうか調べる
//LBRep以外はかならずtrueが返却される
bool MGCurve::cn_continuity(size_t n)const{
	if(type() != MGCURVE_SPLINE)return true;	//LBRep以外に折れはない
	const MGKnotVector &knotvec = knot_vector();
	size_t k, km1, index;
	k = order();
	km1 = k - 1;
	if(k < 2)return true;
	if(km1 <= n){	//C2連続の時orderが4以上かどうか
		//１スパンの時はtrue,それ以上の時はfalseを返却する
		if(k==bdim()){return true;}else{return false;}
	}
	if(!knotvec.locate_multi(km1, k - n, index))return true;
	return false;	//Cn不連続
}

//Construct new curve object by copying to newed area,
//and limitting the parameter range to prange.
//Returned is newed object and must be deleted.
MGCurve* MGCurve::copy_limitted(const MGInterval& prange) const{
	MGCurve* crv=clone();
	crv->limit(prange);
	return crv;
}

//Divide this curve at the designated knot multiplicity point.
int MGCurve::divide_multi(
	MGPvector<MGCurve>& crv_list,	//divided curves will be output.
	int multiplicity	//designates the multiplicity of the knot to divide at.
						//When multiplicity<=0, order()-1 is assumed.
						//When multiplicity>=order(), order() is assumed.
)const{
	crv_list.push_back(this->clone());
	return 1;
}

MGCurve& MGCurve::limit(double t0, double t1){
	MGInterval rng(t0,t1);
	return limit(rng);
}

// 点が曲線上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the curve,
// and 0 if the point is not on the curve.
bool MGPosition::on(
	const MGCurve& curve,	// Curve
	double& t	// Parameter value of the nearest point on the curve.
) const{
	return curve.on(*this, t);
}

// Return curve's parameter value of this point.
// If this point is not on the curve, return the nearest point's parameter
// value on the curve.
double MGPosition::param(const MGCurve& crv) const{
	return crv.param(*this);
}

// Return starting or ending parameter value that is nearer to the param t.
double MGCurve::param_se(double t) const{
	double ts=param_s(),te=param_e();
	return (t-ts)<(te-t) ? ts : te;
}

//Return perpendicular point from a point P,
//given guess starting paramter values.
int MGCurve::perp_guess(
	double t0, double t1,	//parameter range of this.
	const MGPosition& P,	//Point(指定点)
	double tg,				//Guess parameter values of the two curves
	double& t				//Output parameter
)const{

//****Method****
// Let f(t) is this curve, then h=(f-P)**2 should be minimized. 
//If f(t) is the points where the shortest(or longest)
//distance between f and P, then the following conditions are satisfied:
//  dh/dt=2*ft*(f-P)=0.

// Starting with guess parameter t, dt is obtained by solving
// the equation E+dE(t)=0, where E(t)=(df/dt)*(f-P).
//	Then t+dt is the next guess parameter value.
//  dE= (ftt*(f-P)+ft*ft)*dt=-E, where ft is once and 
//  ftt is twice differential of f.
//  If we set A=ftt*(f-P)+ft*ft, dt=-E/A .
//

	bool ranged=t0<t1;
	if(!ranged){
		t0=param_s(); t1=param_e();
	}
	if(tg<t0) tg=t0;
	if(tg>t1) tg=t1;
	double error_sqr=(MGTolerance::wc_zero_sqr())*.25;
	//.25 is multiplied to enforce more strictness of line intersection.
	//Tolerance is made half(.25=.5*.5).

	MGPosition f;
	MGVector ft,ftt,fmP, df;
	double A,E;
	int loop=0, tlow=0, thigh=0;

	t=tg;
	double dt,told=t,tsave;
	while(loop++<16 && tlow<5 && thigh<5){
		eval_all(t,f,ft,ftt);
		fmP=f-P;
		E=ft%fmP;
		A=ftt%fmP+ft%ft;
		if(MGMZero(A)) return 0;

		dt=-E/A; t+=dt;	// Update t.
		df=ft*dt;
		if(df%df<error_sqr){
			int found=in_range(t); t=range(t);
			if(found && ranged){
				double terror=(t1-t0)*MGTolerance::rc_zero();
				found=(t0-terror<=t && t<=t1+terror);
			}
			return found;
		}
		
		tsave=t;
//		if(ranged){
			if(t<t0){
				if(told<t0 && t<=told){ t=t0;return 0;} //If no convergence, return.
				tlow+=1; thigh=0; t=t0;
			}else if(t>t1){
				if(told>t1 && t>=told){ t=t1;return 0;} //If no convergence, return.
				tlow=0; thigh+=1; t=t1;
			}else{tlow=0; thigh=0;}
//		}
		told=tsave;
	}
	return 0;
}

//Return perpendicular points of two curves,
//given guess starting paramter values.
int MGCurve::perp_guess(
	double s0, double s1,		//parameter range of this.
	const MGCurve& curve2,		//2nd curve.
	double t0, double t1,		//parameter range of curve2.
	double sg, double tg,		//Guess parameter values of the two curves
	MGPosition& st				//Output parameter pair of (s,t)
)const{
	const MGCompositeCurve* compoc=dynamic_cast<const MGCompositeCurve*>(&curve2);
	if(compoc){
		if(compoc->perp_guess(t0,t1,*this,s0,s1,tg,sg,st)){
			st.swap(0,1);
			return 1;
		}else
			return 0;
	}

//****Method****
// Let f(s) and g(t) are two curves, then h=(f-g)**2 should be minimized. 
//(f is this curve, and g is curve2.)
//If f(s) and g(t) are the points where the shortest(or longest)
//distance between f and g, then the following conditions are satisfied:
//  dh/ds=2*fs%(f-g)=0.    dh/dt=2*gt%(g-f)=0.

// Starting with guess parameter (s,t), ds and dt are
//  obtained by solving the two equations
//	E+dE(s,t)=0 and F+dF(s,t)=0,
//  where E(s,t)=(df/ds)%(f-g) and F(s,t)=(dg/dt)%(g-f).
//	Then (s+ds, t+dt) is the next guess parameter value.
//	Here dE(s,t) and dF(s,t) are total differentials of E and F.
//  dE= (fss%(f-g)+fs%fs)*ds-fs%gt*dt=-E
//  dF=-fs%gt*ds+(gt%gt+(g-f)%gtt)*dt=-F
//  (fs and gt are once differential of f and g.
//  fss and gtt are twice differential of f and g.)
//  If we set A=fss%(f-g)+fs%fs, B=-fs%gt, and C=gt%gt+(g-f)%gtt,
//  ds=(B*F-C*E)/(A*C-B*B)	dt=(B*E-A*F)/(A*C-B*B) .
//
	bool sranged=s0<s1, tranged=t0<t1;
	if(!sranged){
		s0=param_s(); s1=param_e();
	}
	if(!tranged){
		t0=curve2.param_s(); t1=curve2.param_e();
	}
	if(sg<s0) sg=s0; if(sg>s1) sg=s1;
	if(tg<t0) tg=t0; if(tg>t1) tg=t1;
	double error_sqr=(MGTolerance::wc_zero_sqr())*.25;
	//.25 is multiplied to enforce more strictness of line intersection.
	//Tolerance is made half(.25=.5*.5).

	MGPosition f,g;
	MGVector fs,fss,gt,gtt,fmg, df,dg;
	double A,B,C,E,F,ACmB2;
	int loop=0, slow=0, shigh=0, tlow=0, thigh=0;

	st=MGPosition(sg,tg); double s=sg, t=tg;
	double ds,dt,sold=s,told=t,ssave,tsave;
	while(loop++<16 && slow<5 && shigh<5 && tlow<5 && thigh<5){
		eval_all(s,f,fs,fss); curve2.eval_all(t,g,gt,gtt);
		fmg=f-g;
		E=fs%fmg; F=-(gt%fmg);
		A=fss%fmg+fs%fs; B=-(fs%gt); C=gt%gt-fmg%gtt;
		ACmB2=A*C-B*B;
		if(MGMZero(ACmB2)){
			if(MGMZero(B)) return 0;
			ds=dt=-E*.5;
			dt/=B;
			ds/=A;			
		}else{
			ds=(B*F-C*E)/ACmB2; dt=(B*E-A*F)/ACmB2;
		}
		df=fs*ds; dg=gt*dt;
		// Update s,t.
		s+=ds; t+=dt;	
		if(df%df<error_sqr && dg%dg<error_sqr){
			st(0)=s; st(1)=t;
			int found=in_range(s) && curve2.in_range(t);
			s=range(s); t=curve2.range(t);
			if(found){
				double error_rel=MGTolerance::rc_zero();
				if(sranged){
					double serror=(s1-s0)*error_rel;
					found=(s0-serror<=s && s<=s1+serror);
				}
				if(found){
					if(tranged){
						double terror=(t1-t0)*error_rel;
						found=(t0-terror<=t && t<=t1+terror);
					}
				}
			}
			return found;
		}

		ssave=s; tsave=t;
		//if(sranged){
			if(s<s0){
				if(sold<s0 && s<=sold)
					return 0; //If no convergence, return.
				slow+=1; shigh=0; s=s0;
			}else if(s>s1){
				if(sold>s1 && s>=sold)
					return 0; //If no convergence, return.
				slow=0; shigh+=1; s=s1;
			}else{slow=0; shigh=0;}
		//}
		//if(tranged){
			if(t<t0){
				if(told<t0 && t<=told)
					return 0; //If no convergence, return.
				tlow+=1; thigh=0; t=t0;
			}else if(t>t1){
				if(told>t1 && t>=told)
					return 0; //If no convergence, return.
				tlow=0; thigh+=1; t=t1;
			}else {tlow=0; thigh=0;}
		//}
		sold=ssave; told=tsave;
	}
	return 0;
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//**** NOTE 1 ****
//perpendiculars is the general function of perps, used in perps.
//General users should use function perps, not perpendiculars, since
//perps is optimized for each curve type.
//**** NOTE 2 ****
//perpendiculars can not be used for infinite parameter range curve.
//param_s() and param_e() of both curves must return their finite
//parameter range.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGCurve::perpendiculars(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list;

	double t1s_init=param_s() ,t1s;
	double t1e_init=param_e(), t1e;
	size_t ndiv1=intersect_dnum(); size_t ndiv1m1=ndiv1-1;
	double span1=(t1e_init-t1s_init)/double(ndiv1);

	double t2s_init=crv2.param_s() ,t2s;
	double t2e_init=crv2.param_e(), t2e;
	size_t ndiv2=crv2.intersect_dnum(); size_t ndiv2m1=ndiv2-1;
	double span2=(t2e_init-t2s_init)/double(ndiv2);

	MGPosition ppair;
	//Compute by subdividing both curves to parts.
	t1s=t1s_init;
	for(size_t i=0; i<ndiv1; i++){
		if(i<ndiv1m1) t1e=t1s_init+span1*double(i+1);
		else t1e=t1e_init;
		//subdivide this curve.
		double guess1=(t1s+t1e)/2.;
		t2s=t2s_init;
		for(size_t j=0; j<ndiv2; j++){
			if(j<ndiv2m1) t2e=t2s_init+span2*double(j+1);
			else t2e=t2e_init;
			//subdivide crv2.
			double guess2=(t2s+t2e)/2.;
			if(perp_guess(t1s,t1e,crv2,t2s,t2e,guess1,guess2,ppair))
				list.append(*this,crv2,ppair);
			t2s=t2e;
		}
		t1s=t1e;
	}
	return list;
}

//Compute all foot points of the perpendicular line from point to
//the curve.
// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
// すべて求める。
MGCParam_list MGCurve::perps(
	const MGPosition& P		//Point(指定点)
)const{
	MGCParam_list list(this);

	double t1s_init=param_s() ,t1s;
	double t1e_init=param_e(), t1e;
	size_t ndiv1=intersect_dnum(); size_t ndiv1m1=ndiv1-1;
	double span1=(t1e_init-t1s_init)/double(ndiv1);

	double t;
	//Compute by subdividing this curve to parts.
	t1s=t1s_init;
	for(size_t i=0; i<ndiv1; i++){
		if(i<ndiv1m1) t1e=t1s_init+span1*double(i+1);
		else t1e=t1e_init;
		//subdivide this curve.
		double guess1=(t1s+t1e)/2.;
		if(perp_guess(t1s,t1e,P,guess1,t)) list.append(t);
		t1s=t1e;
	}
	return list;	
}

//Compute all foot points of the perpendicular line from this point to
//a curve.
// ポイントから与曲線へ下ろした垂線の足の，曲線のパラメータ値を
// すべて求める。
MGCParam_list MGPosition::perps(
	const MGCurve& crv2		//Curve
)const{
	return crv2.perps(*this);
}

//Perpendicular points with C1 conitnuity LBRep.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGCurve::perps_withC1LB(
   const MGLBRep& lbC1
)const{
	return perpendiculars(lbC1);
}

//Perpendicular points with SurfCurve
//whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list MGCurve::perps_with_noCompoSC(const MGSurfCurve& curve2)const{
	return perpendiculars(curve2);
}

//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGCurve::perpsSl(
	const MGStraight& sl	//The second curve
)const{
	MGPosition_list list;

	double t1s_init=param_s() ,t1s;
	double t1e_init=param_e(), t1e;
	double range1=t1e_init-t1s_init;
	size_t ndiv1=intersect_dnum(); size_t ndiv1m1=ndiv1-1;
	double dndiv1=double(ndiv1);

	MGPosition ppair;
	//Compute by subdividing this LBRep to parts.
	t1s=t1s_init;
	for(size_t i=0; i<ndiv1; i++){
		if(i<ndiv1m1) t1e=t1s_init+range1*double(i+1)/dndiv1;
		else t1e=t1e_init;
		//subdivide this curve.
		double guess1=(t1s+t1e)/2.;
		if(perp_guess(t1s,t1e,sl,1.,0.,guess1,0.,ppair))
			list.append(*this,sl,ppair);
		t1s=t1e;
	}
	return list;	
}

// Curve上の与ポイントでのCurveの曲率を返却する。
double MGCurve::curvature(double d)const{
	d = range(d);
	double cur;
	MGPosition p; MGVector v1,v2;

	eval_all(d,p,v1,v2);
	double v1_len = v1.len();
	if(MGMZero(v1_len)) return 0.;
	double v1_len3 = v1_len * v1_len * v1_len;
	if(p.sdim()==2) cur = (v1.ref(0)*v2.ref(1) - v1.ref(1)*v2.ref(0))/v1_len3;
	//When sdim==2, curvature has sign, and sdim()=3 curvature is always plus.
	else            cur = (v1*v2).len()/v1_len3;
	return cur;
}

// Curve上の与ポイントでの単位接ベクトルを返却する。
MGUnit_vector MGCurve::direction(double d) const{
	return MGUnit_vector (eval_deriv(d));
}

//Compute direction unit vector of the geometry.
MGUnit_vector MGCurve::direction(const MGPosition& param) const{
	return direction(param(0));
}

// 終点を返却する。
MGPosition MGCurve::end_point() const{
	return eval(param_e());
}

//Obtain an extrapolated PP-Rep curve by the parameter value.
void MGCurve::extrapolated_pp(
	double tau,		//The parameter value at the end of extended point.
					//When tau<param_s(), extension will be done at the starting point.
					//When tau>param_e(), extension will be done at the end point.
	double dk,     //Coefficient of how curvature should vary at the connecting point.
	MGPPRep& pp
)const{
	double t,dtau,t0=param_s(), t1=param_e();
	bool at_start;
	if(tau<t0){
		t=t0; at_start=true; dtau=t0-tau;
	}else{
		t=t1; at_start=false; dtau=tau-t1;
	}
	size_t k=order();
	if(4<k)
		k=4;
	MGPPRep pp2(k,2,sdim());
	pp2.break_point(0)=0.; pp2.break_point(1)=dtau;

	MGPosition P;		//Position data at the extension point.
	MGVector v1,v2,v3;	//i-th derivative data at the extension point.
	eval_all(t,P,v1,v2); v3=eval(t,3);
	double v1len=v1.len();
	if(MGMZero(v1len)){
		//The case when the length of 1st is zero.
		pp2.store_at(0,0,P);
	}else{
		if(at_start)
			v1*=-1.;
		MGUnit_vector g1(v1);
		MGVector v12=v1*v2;
		double v12len=v12.len();
		MGVector g2(0.,0.,0.);
		double curvature,torsion;
		double v1len2=v1len*v1len; double v1len3=v1len2*v1len;
		if(MGMZero(v12len)){
			curvature=torsion=0.;
		}else{
			g2=MGUnit_vector(v12)*g1;
			curvature=v12len/v1len/v1len2;
			torsion=MGDeterminant(v1,v2,v3)/v12len/v12len;
		}
		MGVector g3=g1*g2;//(g1,g2,g3) make Frenet frame.
		double dkdt=-dk*curvature/(v1len*dtau);
			//dkdt=d(curvatue)/dt;how the curvature changes.

		//The following are the Bouquet's formula. See Frene-frame explanation.
		pp2.store_at(0,0,P);
		pp2.store_at(1,0,v1);
		if(k>=3)
			pp2.store_at(2,0,v2);
		if(k>=4)
			pp2.store_at(3,0,(-curvature*curvature*g1+dkdt*g2+curvature*torsion*g3)*v1len3);
	}
	pp=pp2;
}

// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
void MGCurve::eval_all (
		double t,			// パラメータ値
		MGPosition& p,		// 位置
		MGVector& v1,		// 一次微分値
		MGVector& v2		// 二次微分値
		) const{
	t=range(t); p=eval(t,0); v1=eval(t,1); v2=eval(t,2);
}

// 曲線上の与えられたパラメータ値における一次微分値をかえす。
MGVector MGCurve::eval_deriv(double t) const{return eval(t,1);}
									
//Evaluate line data at data point tau.
void MGCurve::eval_line(
	const MGNDDArray& tau,	//Data points.
	MGBPointSeq& value		//Values evaluated. value(i,.)=eval(tau[i]);
)const{
	size_t n=tau.length();
	value.resize(n,sdim());
	for(size_t i=0; i<n; i++) value.store_at(i,eval(tau[i]));
}

// 与えられたパラメータ値に相当する自身上の点を返す。
MGPosition MGCurve::eval_position(double t) const {return eval(t);}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGCurve::evaluate(
	const MGPosition& t,	// Parameter value.
				//t's space dimension is geometry's manifold dimension.
	const size_t* nderiv	//Order of derivative of i-th parameter
				//in nderiv[i].
				//When nderiv=null, nderiv[i]=0 is assumed for all i.
)const{
	if(nderiv) return eval(t(0),nderiv[0]);
	else       return eval(t(0));
}

//Compute Frenet_frame, curvature and torsion in 3D space.
void MGCurve::Frenet_frame(
	double t,			//Input parameter value(パラメータ値)
	MGUnit_vector& T,	//Tangent
	MGUnit_vector& N,	//Principal Normal
	MGUnit_vector& B,	//Binormal
	double& curvature,	//Curvature is always >=0.
	double& torsion
)const{
	MGVector T2,N2,B2;
	MGVector v2;
	Frenet_frame2(t,v2,T2,N2,B2);
	double mzero=MGTolerance::mach_zero();
	double verocity=T2.len();
	if(verocity<=mzero){
		curvature=torsion=0.;
		T=MGVector(1.,0.,0.); N=MGVector(0.,1.,0.); B=MGVector(0.,0.,1.);  
	}else{
		T=T2;
		double B2_len=B2.len();
		if(B2_len<=mzero){
			curvature=torsion=0.;
			T.orthonormal(v2,N,B);
		}else{
			B=B2; N=N2;
			curvature=B2_len/(verocity*verocity*verocity);
			MGVector v3=eval(t,3);
			torsion=MGDeterminant(T2,v2,v3)/(B2_len*B2_len);
		}
	}
}

//Compute Frenet_frame, curvature and torsion in 3D space.
void MGCurve::Frenet_frame2(
	double t,			//Input parameter value(パラメータ値)
	MGVector& V2,//2nd derivative at t.
	MGVector& T,//Tangent
	MGVector& N,//Principal Normal
	MGVector& B	//Binormal
)const{
	MGPosition P;
	eval_all(t,P,T,V2);
	B=T*V2;
	N=B*T;
}

//Test if this curve has the same direction with curve2 at the point s(of this)
// and t(of curve2).
//Function's return value is true if they have the same direction.
//"same direction" means their tangent vectors have the angle less than 90 degree.
bool MGCurve::has_same_direction_at(
	double s,
	const MGCurve& curve2,
	double t
)const{
	MGVector tan1=eval(s,1);
	MGVector tan2=curve2.eval(t,1);
	return tan1%tan2>0.;
}

//Test if input parameter value is inside parameter range of the line.
bool MGCurve::in_range(double t)const{
	const double t1=param_s(), t2=param_e();
	double error=(t2-t1)*MGTolerance::rc_zero();
	return (t>=t1-error && t<=t2+error);
}

//Test if input parameter value is inside parameter range of the line.
bool MGCurve::in_range(const MGPosition& param) const{
	return in_range(param[0]);
}

//Test if this cure is co-planar with the 2nd curve curve2.
//MGPlane expression will be out to plane if this is co-planar.
//Function's return value is true if co-planar.
bool MGCurve::is_coplanar(const MGCurve& curve2, MGPlane& plane)const{
	if(!is_planar(plane))
		return false;
	MGPlane plane2;
	if(!curve2.is_planar(plane2))
		return false;

	return plane==plane2;
}

//Test if the input parameter t is the start point parameter or not.
bool MGCurve::is_startpoint_parameter(double t)const{
	return MGREqual_base(param_s(),t,param_span());
}

//Test if the input parameter t is the start point parameter or not.
bool MGCurve::is_endpoint_parameter(double t)const{
	return MGREqual_base(param_e(),t,param_span());
}

//Test if the vector from P to this->eval(t) is perpendicular to
//the tangent of this curve at t.
bool MGCurve::is_perpendicular(const MGPosition& P, double t)const{
	MGVector V=P-eval(t);
	if(MGAZero(V.len()))
		return true;

	return V.orthogonal(eval(t,1));
}

//Test if this cure is linear or not, that is, is straight or not.
//MGStraight expression will be out to straight if this is linear or not.
//Function's return value is true if linear.
bool MGCurve::is_linear(MGStraight& straight)const{
	size_t ndiv=intersect_dnum();
	double ts=param_s(), te=param_e();
	double tmid=(ts+te)*.5;
	MGUnit_vector dir=eval(tmid,1);
	MGPosition Ps=eval(tmid);
	straight=MGStraight(dir,Ps);
	double dt=(te-ts)/double(ndiv), para;
	for(size_t i=0; i<=ndiv; i++){
		double t=te-dt*double(i);
		MGPosition P=eval(t);
		if(!straight.on(P,para))
			return false;
	}
	return true;
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGCurve::is_planar(MGPlane& plane)const{
	MGLBRep acrv(*this);
	return acrv.is_planar(plane);
}

// Curve と Curve の交点を求める。
//***Caution***intersect can be used only for finite curve, i.e.
//parameter range of the computation is only from param_s() to param_e().
//For example, "intersect" cannot be applied to infinite straight line.
MGCCisect_list MGCurve::intersect(const MGCurve& l2) const{
	MGCCisect_list list(this, &l2);
	if(!has_common(l2))
		return list;

	double error=MGTolerance::wc_zero();
	double tlen1, tlen2, tm1, tm2;
	size_t len_stack=200, max_add=5, current_add=1;;
	std::vector<double> stack(len_stack*4);
//Actually, stack is used as: stack[i][4], where
// stack[i][0] - stack[i][1] is the own line parameter range(min and max).
// stack[i][2] - stack[i][3] is l2 line parameter range.

	double ts1=param_s(), te1=param_e(); 
	double tmax1=(te1-ts1)/intersect_dnum();	//Maximum parameter length
			//curve will be subdivided until the parameter length<=tmin1.
	double ts2=l2.param_s(), te2=l2.param_e();
	double tmax2=(te2-ts2)/l2.intersect_dnum();	//Same as tmax1
	//[ts1,te1] and [ts2,te2] are current spans to check.

	stack[0]=ts1; stack[1]=te1;
	stack[2]=ts2; stack[3]=te2;
	size_t stack_id=1;	// stack_id indicates next available id for stack, i.e.
						// how many data in stack.
 
	MGPosition t12; MGBox box1,box2;
	while(stack_id){
		stack_id -=1;	// Pop the stack data.
		ts1=stack[stack_id*4];		te1=stack[stack_id*4+1];
		ts2=stack[stack_id*4+2];	te2=stack[stack_id*4+3];
		box1=box_limitted(MGInterval(ts1,te1));
		box2=l2.box_limitted(MGInterval(ts2,te2));
loop30:
		//Check of box boundary override-ness.
		if(!((box1&box2).empty())){
			tm1=(ts1+te1)*.5; tm2=(ts2+te2)*.5;
			tlen1=(te1-ts1)/tmax1; tlen2=(te2-ts2)/tmax2;
			if(tlen1<=1. && tlen2<=1.){
			//Compute using perp_guess
				if(perp_guess(ts1,te1,l2,ts2,te2,tm1,tm2,t12)){
					tm1=t12.ref(0); tm2=t12.ref(1);
					MGPosition point1=eval(tm1), point2=l2.eval(tm2);
					if((point1-point2).len()<=error)
						list.append(MGCCisect(point1,tm1,tm2));
				}
			}else{
			//Subdivide and put on stack.
				if(stack_id >= len_stack){
					if(current_add>max_add) break;//Halt the computation.
					else{
						len_stack+=200; stack.resize(len_stack*4);
						current_add+=1;
					}
				}
				stack[stack_id*4+1]=te1; stack[stack_id*4+3]=te2;
				if(tlen1>=tlen2){	//Subdivide own curve.
					stack[stack_id*4]=tm1; stack[stack_id*4+2]=ts2;
					te1=tm1;
					box1=box_limitted(MGInterval(ts1,te1));
				}else{				//Subdivide l2.
					stack[stack_id*4]=ts1; stack[stack_id*4+2]=tm2;
					te2=tm2;
					box2=l2.box_limitted(MGInterval(ts2,te2));
				}
				stack_id +=1;
				goto loop30;
			}
		}
	}		
	
	return list;
}

class MGCurve_isec1DDrive{
	const MGCurve* m_curve;	//SurfCurve to evaluate.
	size_t m_cod;				//coordinate kind of m_f.
	double m_f;					//coordinate value of intersect_1D.
public:
	MGCurve_isec1DDrive(const MGCurve* scrv, size_t cod, double f)
		:m_curve(scrv),m_cod(cod), m_f(f){;};
	double operator()(double t){
		MGVector P=m_curve->eval(t);
		return P[m_cod]-m_f;
	}
};

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
//isect_1Dから座標種 coordinateに平行な平面との交点を求める。
MGCParam_list MGCurve::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
)const{
	size_t ndiv=intersect_dnum();
	double t0=param_s(), t1=param_e();
	double delta=(t1-t0)/double(ndiv);
	double error=MGTolerance::wc_zero();

	//construct function object to evaluate the distance between
	//f and eval(t)[coordinate].
	MGCurve_isec1DDrive diff(this,coordinate,f);

	double tpre,taft=t0;
	double dpre, daft=diff(taft);
	MGCParam_list list(this);

	//Iterate by checking singned distance from the coordinate vaue to the value f.
	//When dpre and daft have different signs, an intersection point must 
	//lie between tpre and taft.
	size_t i=0;
	while(i<=ndiv){
		dpre=daft; tpre=taft;
		while(fabs(dpre)<=error){
			list.append(tpre);
			if(i>=ndiv) break;
			i++; tpre=t0+delta*double(i);
			dpre=diff(tpre);
		}
		if(i>=ndiv) break;
		i++; taft=t0+delta*double(i);
		daft=diff(taft);
		if(fabs(daft)<=error) continue;
		else if(dpre*daft<0.){
		//Now there exists a solution between tpre and taft.
			int ier;
			double x=mgNlbit(diff, tpre,taft, error, 20, ier);
			list.append(x);
		}
	}

	return list;
}

// Curve と Curve の交点を求める。
MGCCisect_list MGCurve::intersect_brute_force(const MGCurve& l2) const{
	size_t i; double tlen1, tlen2, tm1, tm2;
	double errbnd[3]={1., 1.2, 1.35};
	size_t len_stack= 200, max_add=5, current_add=1;;
	std::vector<double> stack(len_stack*4);
//Actually, stack is used as: stack[i][4], where
// stack[i][0] - stack[i][1] is the own line parameter range(min and max).
// stack[i][2] - stack[i][3] is l2 line parameter range.

	size_t ncd=sdim();
	const size_t ncd2=l2.sdim();
	if(ncd2>ncd) ncd=ncd2;

	double error=MGTolerance::wc_zero()*.7;
	size_t errbnd_id=ncd; if(ncd>3) errbnd_id=3;
	double error2=error*errbnd[errbnd_id-1];

	double ts1=param_s(), te1=param_e(); 
	double ts2=l2.param_s(), te2=l2.param_e();
	double xs1=ts1-1., xs2=ts2-1.;
	//[ts1,te1] and [ts2.te2] are current spans to check.
	//xs1 and xs2 are previous solution obtained.

	MGBox box1, box2; MGPosition p1low, p1high, p2low, p2high;

	MGCCisect_list list(this, &l2);

	stack[0]=ts1; stack[1]=te1;
	stack[2]=ts2; stack[3]=te2;
	size_t stack_id=1;	// stack_id indicates next available id for stack, i.e.
						// how many data in stack.
 
	while(stack_id){
		stack_id-=1;	// Pop the stack data.
		ts1=stack[stack_id*4];		te1=stack[stack_id*4+1];
		ts2=stack[stack_id*4+2];	te2=stack[stack_id*4+3];
		box1=box_limitted(MGInterval(ts1,te1));
		box2=l2.box_limitted(MGInterval(ts2,te2));
loop30:
		//Check of box boundary override-ness.
		if((box1&box2).empty())	goto loop_end;
		p1low=box1.low(), p1high=box1.high();
		p2low=box2.low(), p2high=box2.high();
		tlen1= tlen2 =0.;
//		double f1l,f1h, f2l,f2h;////////////////////
		for(i=0; i<ncd; i++){
//			 f1h=p1high.ref(i); f1l=p1low.ref(i);	//////////
//			 f2h=p2high.ref(i); f2l=p2low.ref(i);////////////
			tlen1+=p1high.ref(i)-p1low.ref(i);
			tlen2+=p2high.ref(i)-p2low.ref(i);
		}					//tlen1,2 are sum of each axis length of the box.
		tm1=(ts1+te1)*.5; tm2=(ts2+te2)*.5;
		if(tlen1<=error2 && tlen2<=error2){
			if(fabs(tm1-xs1)>error && fabs(tm2-xs2)>error){
				MGPosition point=eval(tm1);
				list.append(MGCCisect(point,tm1,tm2));
				xs1=tm1; xs2=tm2;
			}
		}
		else{
			if(stack_id >= len_stack){
				if(current_add>max_add) break;
				else{
					len_stack+=200; stack.resize(len_stack*4);
					current_add+=1;
				}
			}
			stack[stack_id*4+1]=te1; stack[stack_id*4+3]=te2;
			if(tlen1>=tlen2){	//Subdivide own curve.
				stack[stack_id*4]=tm1; stack[stack_id*4+2]=ts2;
				if(tm1<xs1 && xs1<=te1 && ts2<xs2 && xs2<=te2) goto loop_end;
				te1=tm1;
				box1=box_limitted(MGInterval(ts1,te1));
			}
			else{				//Sbudivide l2 curve.
				stack[stack_id*4]=ts1; stack[stack_id*4+2]=tm2;
				if(ts1<xs1 && xs1<=te1 && tm2<xs2 && xs2<=te2) goto loop_end;
				te2=tm2;
				box2=l2.box_limitted(MGInterval(ts2,te2));
			}
			stack_id +=1;
			goto loop30;
		}
loop_end:;
	}
	
	return list;
}

// MGRLBRep と Curve の交点を求める。
MGCCisect_list MGCurve::isect(const MGRLBRep& curve2)const{
	return intersect(curve2);
}

//Intersection of MGCurve and MGEllipse.
MGCCisect_list MGCurve::isect(const MGEllipse& curve2)const{
	return intersect(curve2);
}

//Intersection point of curve and MGLBRep.
MGCCisect_list MGCurve::isect(const MGLBRep& curve2) const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection of MGCurve and MGBSumCurve.
MGCCisect_list MGCurve::isect(const MGBSumCurve& curve2)const{
	return intersect(curve2);
}

MGCCisect_list MGCurve::isect(const MGTrimmedCurve& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

MGCCisect_list MGCurve::isect(const MGCompositeCurve& curve2)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

MGCSisect_list MGCurve::isect(const MGFace& f)const{
	MGCSisect_list list;
	const MGSurface* srf=f.surface();
	if(!srf)
		return list;
	if(!f.has_common(*this))
		return list;

	list=isect(*srf);
	MGCSisect_list::CSiterator	i=list.begin(), iend=list.end(), i1;
	while(i!=iend){
		i1=i; i1++;
		if(!f.in_range((*i).param_surface()))
			list.removeAt(i);
		i=i1;
	}
	return list;
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list MGCurve::isect_withC1LB(const MGLBRep& curve2)const{
	return intersect(curve2);
}

//isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list MGCurve::isect_with_noCompoSC(const MGSurfCurve& curve2)const{
	return intersect(curve2);
}

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
MGCParam_list MGCurve::isect_1D(
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
) const{
	const MGInterval& rng=box()[coordinate];
	double error=MGTolerance::wc_zero();
	if(rng[0]>f+error || rng[1]<f-error) return MGCParam_list();
	return intersect_1D(f,coordinate);
}

//Intersection of a curve and a plane
MGCSisect_list MGCurve::intersect_with_plane(const MGPlane& surf) const{
	double g[4]={surf.normal().ref(0),	surf.normal().ref(1) 
				,surf.normal().ref(2), surf.distance()};
	std::auto_ptr<MGCurve> crv1D=oneD(g);// std::cout<<(*crv1D)<<std::endl;//////////

	MGCSisect_list list(this,&surf);
	double error=MGTolerance::wc_zero();
	const MGInterval& minmax=crv1D->box()[0];
	if(minmax[0]>error || minmax[1]<-error)
		return list;

	MGCParam_list clist=crv1D->isect_1D(0.);
	size_t inum=clist.entries();
	MGPosition p,uvpl;
	for(size_t i=0; i<inum; i++){
		double tt=clist.removeFirst();
		p=eval(tt);
		uvpl=surf.uv(p);
		const MGCSisect isect(surf.eval(uvpl),tt,uvpl);
		list.append(isect);
	}

	return list;
}

//Transform the coordinates of boundary of this geometry so that
//new coordinate of boundary is the same coordinate as the new one of
//this geometry after negate() of this geometry is done.
//That is, boundary coordinates are of parameter world of this geometry.
void MGCurve::negate_transform(MGGeometry& boundary) const{
	assert(dynamic_cast<MGPoint*>(&boundary));
	MGPoint* param=dynamic_cast<MGPoint*>(&boundary);
	double t=param->position().ref(0);
	t=negate_param(t);
	boundary=MGPoint(MGPosition(1,&t));
}

// 点がCurve上にあるか調べる。Curve上であれば，そのパラメータ値を，
// そうでなくても最近傍点のパラメータ値を返す。
bool MGCurve::on(
	const MGPosition& point,	// 指定点
	double& t				// Parameter value will be returned.
)const{
	size_t i; MGCParam_list tlist;
	double len;
	double tolsqr=MGTolerance::wc_zero_sqr();

	//1.Test about end points.
	double t0=param_s(),t1=param_e();
	MGVector dif0=point-eval(t0), dif1=point-eval(t1);
	double len0=dif0%dif0, len1=dif1%dif1;
	if(len0<len1){
		t=t0; len=len0;
	}else{
		t=t1; len=len1;
	}
	if(len<=tolsqr)
		return true;

	t0=(t0+t1)*0.5;
	dif0=point-eval(t0); len0=dif0%dif0;
	if(len0<len){t=t0; len=len0;}
	if(len<=tolsqr)
		return true;

	//2.Compute the point on the line by obtaining intersection points
	//   with each axis.
	double ti,leni;
	for(i=0; i<sdim(); i++){
		tlist=isect_1D(point.ref(i),i);
		size_t inum=tlist.entries();
		for(size_t j=0; j<inum; j++){
			ti=tlist.removeFirst();
			MGVector dif(point-eval(ti));
			leni=dif%dif;
			if(leni<len){t=ti; len=leni;};
		}
	}
	if(len<=tolsqr)
		return true;

	//3. Compute nearest point.
	if(perp_guess(1.,0.,point,t, ti)){
		MGVector dif(point-eval(ti));
		leni=dif%dif;
		if(leni<len){t=ti; len=leni;}
	}
	
	return len<=tolsqr;
}

//Test if given point is on the geometry or not. If yes, return parameter
//value of the geometry. Even if not, return nearest point's parameter.
// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the geometry,
// and 0 if the point is not on the geometry.
bool MGCurve::on(
	const MGPosition& P,//Point(指定点)
	MGPosition&	param	//Parameter of the geometry(パラメータ)
)const{
	double t;
	bool is_on=on(P,t);
	param=MGPosition(1,&t);
	return is_on;
}

// Return the curve's parameter value of the input point p.
// If input point is not on the curve, return the nearest point's parameter
// value on the curve.
double MGCurve::param(const MGPosition& p)const{
	double t;
	on(p, t);
	return t;
}

//Obtain parameter space error.
double MGCurve::param_error()const{
	return MGTolerance::rc_zero()*param_span();
}

//Return parameter range of the curve(パラメータ範囲を返す)
MGInterval MGCurve::param_range()const{
	return MGInterval(param_s(), param_e());
}

//Return parameter range of the geometry(パラメータ範囲を返す)
MGBox MGCurve::parameter_range()const{
	MGInterval itrvl=param_range();
	return MGBox(1,&itrvl);
}

MGPosition_list MGCurve::perps(
	const MGRLBRep& crv2	//The second curve
)const{
	return perpendiculars(crv2);
}

MGPosition_list MGCurve::perps(
	const MGEllipse& crv2	//The second curve
)const{
	return perpendiculars(crv2);
}

//Perpendicular points with C1 conitnuity LBRep.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGCurve::perps(
	const MGLBRep& crv2		//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGCurve::perps(
	const MGSurfCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGCurve::perps(
	const MGBSumCurve& crv2	//The second curve
)const{
	return perpendiculars(crv2);
}

MGPosition_list MGCurve::perps(
	const MGCompositeCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGCurve::perps(
	const MGTrimmedCurve& crv2	//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

//Compute a foot point of the perpendicular line from point p to
//the curve. If more than one points are found, return nearest one.
// 指定点からの自身への垂線の足とパラメータ値を返す。
// Function's return value is if point is obtained(1) or not(0)
int MGCurve::perp_point(
	const MGPosition &point,	// 与ポイント
	double &t,					// パラメータ値
	const double *g				// guess parameter value.
)const{
	MGCParam_list tlist=perps(point);
	size_t nump=tlist.entries();
	if(nump){
		MGCParam_list::Citerator i,j;
		i=j=tlist.begin();
		double len, lent;
		if(nump>1){
			if(g){	//Compute nearest param to *g
				double param=range(*g);
				len=fabs(param-(*i));
				for(++i; i!=tlist.end(); i++){
					lent=fabs(param-(*i));
					if(lent<len){ j=i; len=lent;};
				}
			}
			else{  //Find minimum length point from p.
				len=(point-eval(*i)).len();
				for(++i; i!=tlist.end(); i++){
					lent=(point-eval(*i)).len();
					if(lent<len){ j=i; len=lent;};
				}
			}
		}
		t=*j;
	}
	return nump;
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGCurve::pick_closest(const MGStraight& sl)const{
	MGPosition param2=closest(sl);
	MGPosition param(1,param2.data());
	return param;
}

//Round t into curve's parameter range.
// 入力パラメータをパラメータ範囲でまるめて返却する。
double MGCurve::range(double t)const{
	double t1=param_s();
	if(t<t1) t=t1;
	else{
		double t2=param_e();
		if(t>t2) t=t2;
	}
	return t;
}

///Approximate this curve by a polyline and output to lb2.
///The tolerance of the approximation is error.
void MGCurve::polygonize(
	double error,	///<tolerance allowed for the approximation
	MGLBRep& lb2	///<Obtained polyline will be output as an MGLBRep of order2.
)const{
	double saved=MGTolerance::set_line_zero(error);
	lb2=MGLBRep(*this,2);
	MGTolerance::set_line_zero(saved);
}

//Round t into geometry's parameter range.
// 入力パラメータをパラメータ範囲でまるめて返却する。
MGPosition MGCurve::range(const MGPosition& prange)const{
	double t=prange(0);
	t=range(t);
	return MGPosition(1,&t);
}

// 指定点をとおり指定方向ベクトルを持つ直線の回りを指定角度
// 回転させて自身の曲線とする。
MGCurve& MGCurve::rotate_self (
	const MGVector & v,
	double d,
	const MGPosition & p
){
	MGTransf t(3); t.set_rotate_3D(v, d, p);// 指定された変換を作成する。
	*this *= t;	// 自身を変換する。
	return *this;
}

//始点を返却する。
MGPosition MGCurve::start_point() const{
	return eval(param_s());
}
