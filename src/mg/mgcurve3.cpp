/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Unit_vector.h"
#include "mg/Transf.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/Gausp.h"
#include "mg/TrimmedCurve.h"
#include "mg/CCisect_list.h"
#include "mg/nlbit.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGCurve
// Implementation of MGCurve.

//The class for function object for mgGausp to compute curvilinear integral
// of the line, compute sqrt(sum(1st deriv of axis i)).
class MGCurve_curvilinear_integral{
	const MGCurve* m_curve;
public:
	MGCurve_curvilinear_integral(const MGCurve* curve):m_curve(curve){;};
	double operator()(double t)const;
};

double MGCurve_curvilinear_integral::operator()(double t)const{
	MGVector pos=m_curve->eval(t,0);
	MGVector deri=m_curve->eval(t,1);
	return (pos.ref(0)*deri.ref(1)-deri.ref(0)*pos.ref(1));
}

double MGCurveLenParamDrive::operator()(double t)const{
	double len=m_curve->length(m_ts, t);
	return len-m_len;
}

//線積分を求める。
//Compute curvilinear integral of the 1st two coordinates.
//This integral can be used to compute area sorrounded by the curve.
double MGCurve::curvilinear_integral(double t1, double t2) const{

	//Get interval.
	double s1=t1, s2=t2;
	if(s1>s2){ s1=t2; s2=t1;}
	s1=range(s1); s2=range(s2);

	//Get knot vector of the curve.
	const MGKnotVector& kv=knot_vector();
	int i1=kv.locate(s1);
	int i2=kv.locate(s2);

	MGCurve_curvilinear_integral clFunc(this);
	double len;
	if(i1==i2) len=mgGausp(clFunc,s1,s2);
	else{
		double len_span,u1,u2;
		u2=kv(i1+1);
		len=mgGausp(clFunc,s1,u2);
		for(int i=i1+1; i<i2; i++){
			u1=kv(i); u2=kv(i+1);
			if(u1<u2){
				len_span=mgGausp(clFunc,u1,u2);
				len +=len_span;
			}
		}
		len_span=mgGausp(clFunc,u2,s2);
		len +=len_span;
	}

	if(t1>t2) len=-len;
	return len;
}

//Compute mean length of 1st derivative vector.
double MGCurve::deriv_length()const{
	double t0=param_s(), t1=param_e();
	double tmid=(t0+t1)*.5;
	double len=eval(t0,1).len();
	len+=eval(t1,1).len();
	len+=eval(tmid,1).len();
	return len/3.;
}

/**
 *  @brief eval_discrete_deviationの下請け関数
 *  @param curve1 face1側エッジのワールドカーブ(must be trimmed into evaluatin range).
 *  @param curve2 face2側エッジのワールドカーブ(must be trimmed into evaluatin range).
 */
void deviation(
	const MGCurve&    curve1,
	const MGCurve&    curve2,
	int npoint,		//num of discrete points.
	std::vector<MGPosition>& sts
){
	MGPosition P1S = curve1.start_point();
	MGPosition P2S = curve2.start_point(), P2E = curve2.end_point();
	bool same_direction=(P1S-P2S).len()<=(P1S-P2E).len();

	double s=curve1.param_s(), ds=curve1.param_span()/npoint;
	double t, dt=curve2.param_span()/npoint;
	if(same_direction){
		t = curve2.param_s();
	}else{
		t = curve2.param_e();
		dt *= -1.;
	}
	sts.push_back(MGPosition(s,t));

	// pos1, pos2の、それぞれの最初と最後の要素以外のすべての要素を計算する
	for(int i = 1; i < npoint; i++){
		s += ds; t += dt;
		MGPosition P1=curve1.eval(s);
		if(!curve2.perp_guess(0.,-1., P1, t, t)){
			t = curve2.closest(P1);
		}
		sts.push_back(MGPosition(s,t));
	}

	// pos1, pos2それぞれの最後の点を決定する
	s=curve1.param_e();
	if(same_direction){
		t=curve2.param_e();
	}else{
		t=curve2.param_s();
	}
	sts.push_back(MGPosition(s,t));
}

//Evaluate deviations of two curves(this and curve2) at npoint discrete points.
//(1)Search the common curve spans which have the distance within tolerance.
//(2)Compute the nearest points from npoint discrete points of this to curve2.
//Let sti=sts[i], then
//sti[0] is this curve's parameter value s, and sti[1] is the parameter value t
//of curve2 which is the nearest point from the point s.
//If this and curve2 have the minimum distance more than tolerance,
//sts.size()==1 and sts[0] is the minimum distance points of this and curve2.
void MGCurve::eval_discrete_deviation(
	const MGCurve& curve2,
	std::vector<MGPosition>& sts,
	int npoint,		//indicates how many discrete points be obtained.
	double tolerance	//tolerance to get two edge to compute deviation.
)const{
	if(tolerance < MGTolerance::line_zero()){
		tolerance= MGTolerance::line_zero();
	}

	MGPosition st=closest(curve2);
	MGPosition P1=eval(st[0]), P2=curve2.eval(st[1]);
	if(P1.distance(P2)>tolerance){
		sts.push_back(st);
		return;
	}

	double tol_saved = MGTolerance::set_line_zero(tolerance);
	std::vector<double> stspans;
	MGCCisect_list isects;
	int retcode = common(curve2, stspans,isects);
	MGTolerance::set_line_zero(tol_saved);

	if(retcode <= 0){//When common failed.
		sts.push_back(st);
		return;
	}else if(retcode==2){//When only intersections are obtained.
		MGCCisect_list::iterator i=isects.begin(), ie=isects.end();
		for(;i!=ie;i++){
			sts.push_back(MGPosition(i->param1(), i->param2()));
		}
		return;
	}
	
	size_t j, nspan=stspans.size()/4;
	double slen=stspans[1]-stspans[0];
	for(j=1; j<nspan; j++){
		slen+=stspans[4*j+1]-stspans[4*j];
	}

	for(j=0; j<nspan; j++){
		size_t j4=j*4;
		double s0 = stspans[j4], s1=stspans[j4+1], t0=stspans[j4+2], t1=stspans[j4+3];
		/*if(s0>s1){
			double save=s0; s0=s1; s1=save;
		}
		if(t0>t1){
			double save=t0; t0=t1; t1=save;
		}*/

		// 共線部分曲線を取り出す
		MGTrimmedCurve crv1trmd(*this,s0, s1), crv2trmd(curve2,t0, t1);
		deviation(crv1trmd,crv2trmd,int(double(npoint)*(s1-s0)/slen),sts);
	}
}

MGKnotVector& MGCurve::knot_vector(){
	const MGCurve* crv=const_cast<const MGCurve*>(this);
	const MGKnotVector& t=crv->knot_vector();
	return *(const_cast<MGKnotVector*>(&t));
}

// 与えられたパラメータ値間の曲線の長さを返す。
// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
//Return line length from t1 to t2.
double MGCurve::length(double t1, double t2) const{
	double s1=t1, s2=t2;
	if(s1>s2){ s1=t2; s2=t1;}
	s1=range(s1); s2=range(s2);

	//Get knot vector of the curve.
	const MGKnotVector& kv=knot_vector();
	int i1=kv.locate(s1);
	int i2=kv.locate(s2);

	double u1,u2;
	MGCurveLengthDrive lenFunc(this);
	double len;
	if(i1==i2) len=mgGausp(lenFunc,s1,s2);
	else{
		double len_span;
		u2=knot(i1+1);
		len=mgGausp(lenFunc,s1,u2);
		for(int i=i1+1; i<i2; i++){
			u1=knot(i); u2=knot(i+1);
			len_span=mgGausp(lenFunc,u1,u2);
			len +=len_span;
		}
		len_span=mgGausp(lenFunc,u2,s2);
		len +=len_span;
	}
	if(t1>t2) len=-len;
	return len;
}

// パラメータで示される点 t_in から指定距離 len はなれた点のパラメータ
// を返す。
double MGCurve::length_param(double t_in, double len) const{
	double ts=range(t_in);
	double tlimit1=param_s();
	double tlimit2=param_e();
	double tlimit = (len>=0.) ? tlimit2:tlimit1;
	MGVector deriv1=eval(ts,1); double dl1=deriv1.len();
	if(MGMZero(dl1)) return tlimit;
	double d=range(ts+len/dl1);	//d is approximate solution.

	// Check if we have solution.
	double len_limit=length(ts,tlimit);
	if(fabs(len)>=fabs(len_limit)) return tlimit; 
									//Solution is outside range;
	double t1,t2,len1,lenm;
	if(ts<=tlimit){
		t1=ts; t2=tlimit;
		len1=-len;
	}else{
		t2=ts; t1=tlimit;
		len1=len_limit-len;
	}
	lenm=length(ts,d)-len;
	if(len1*lenm<=0.) t2=d; else t1=d;

	//Now there exists a solution between t1 and t2.
	int itr=30, ier; 
	double eps=MGTolerance::wc_zero();	
	MGCurveLenParamDrive clen(this,len,ts);
	return mgNlbit(clen, t1,t2, eps, itr, ier);
}

//Compute curvatur in 3D space, ie, the value is not negative.
double MG_Curvature(const MGVector& v1, const MGVector& v2){
	double v1_len = v1.len();
	if(MGMZero(v1_len)) return 0.;
	double v1_len3 = v1_len * v1_len * v1_len;
	return (v1*v2).len()/v1_len3;
}

//Compute torsion.
MGEXTERN double MG_Torsion(
		const MGVector& v1,		//First derivative.
		const MGVector& v2,		//Second derivative.
		const MGVector& v3)		//Third derivative.
{
	MGVector v12=v1*v2;
	double v12_len2=v12%v12;
	if(MGMZero(v12_len2)) return 0.;
	else return MGDeterminant(v1,v2,v3)/v12_len2;
}

double MGCurveLengthDrive::operator()(double t)const{
	MGVector deriv=m_curve->eval_deriv(t);
	return sqrt(deriv%deriv);
}

//Generate arrow data of the tangent at the parameter value t of the curve.
//data[0] is the origin of the arrow, data[1] is the top of the arrow,
//data[2], [3] are two bottoms of arrowhead.
void MGCurve::arrow(double t, MGPosition data[4])const{
	const double arrow_length=.1//of total length of the curve)
				, head_length=.3;//of arrow_length.

	MGVector v1,v2;
	eval_all(t,data[0],v1,v2);
	v1*=param_span()*arrow_length;
	one_arrow(data[0],v1,v1*v2*v1,data[1],data[2],data[3]);
}

//Generate arrow data from (root, vecx, vecy).
void one_arrow(
	const MGPosition& root,	//root of the arrow
	const MGVector& vecx,	//the vector from the root to the head of the arrrow
	const MGUnit_vector& vecy,//vecy that is normal to the vector from root to head
	MGPosition& head,		//head of the arrow will be returned.
	MGPosition& headtail1,	//two tail of arrowhead line segments will be returned.
	MGPosition& headtail2
){
	const double arrow_length=.1//of total length of the curve)
				, head_length=.3;//of arrow_length.

	head=root+vecx;
	MGVector head_rootx=vecx*head_length,
				head_rooty=vecy*(.5*head_length*vecx.len());
	headtail1=head-head_rootx+head_rooty;
	headtail2=head-head_rootx-head_rooty;
}

//Round the parameter t into this parameter range.
double MGCurve::param_round_into_range(double t)const{
	if(t<param_s()) return param_s();
	if(t>param_e()) return param_e();
	return t;
}

class MGtangent_guess_cangle{
	const MGCurve& m_curve;
	const MGPosition& m_P;
public:
	MGtangent_guess_cangle(
		const MGCurve& curve,
		const MGPosition& P
	):m_curve(curve),m_P(P){;};

	double operator()(double t)const{
		MGPosition B;
		MGVector tan, deri2;
		m_curve.eval_all(t,B,tan,deri2);
		return (m_P-B).cangle(tan*deri2*tan);
	}
};

#define MAX_ITR 10
//Return tangent point from a point P,
//given guess starting paramter tg.
//   tangent_guess=true if perpendicular points obtained,
//   tangent_guess=false if perpendicular points not obtained,
int MGCurve::tangent_guess(
	double t0, double t1,	//parameter range of this.
	const MGPosition& P,	//Point(指定点)
	double tg,				//Guess parameter values of the two curves
	double& t				//Output parameter
)const{
	//initial check.
	MGPosition BP;
	MGVector tan, deri2;
	eval_all(tg,BP,tan,deri2);
	double mzero=MGTolerance::mach_zero();
	if(tan.len()<=mzero) return 0;
	if(deri2.len()<=mzero){
		if(tan.is_collinear(P-BP)){
			t=tg;
			return 1;
		}
		return 0;
	}

	if(t0>=t1){
		t0=param_s(); t1=param_e();
	}
	if(tg<t0) tg=t0;
	if(tg>t1) tg=t1;
	double error=MGTolerance::rc_zero();

	int ndiv=intersect_dnum();
	double span=(param_e()-param_s())/double(ndiv);

	int no_converged;
	MGtangent_guess_cangle cang(*this,P);
	double cangg=cang(tg);
	double tn1=tg+span; if(tn1>t1) tn1=t1;
	double cangn1=cang(tn1);
	if(cangg*cangn1<=0.){
		t=mgNlbit(cang,tg,tn1,error,MAX_ITR,no_converged);
		if(no_converged) return 0;
		return 1;
	}

	double tn2=tg-span; if(tn2<t0) tn2=t0;
	double cangn2=cang(tn2);
	if(cangg*cangn2<=0.){
		t=mgNlbit(cang,tg,tn2,error,MAX_ITR,no_converged);
		if(no_converged) return 0;
		return 1;
	}

	if(fabs(cangn1)>fabs(cangn2)) span*=-1.;
	for(int loop=0; loop<ndiv; loop++){
		tn1=tg+span;
		if(tn1>t1){ tn1=t1; loop=ndiv;}
		if(tn1<t0){ tn1=t0; loop=ndiv;}
		cangn1=cang(tn1);

		if(cangg*cangn1<=0.){
			t=mgNlbit(cang,tn1-span,tn1,error,MAX_ITR,no_converged);
			return no_converged-1;
		}
	}
	return 0;
}

//Trim the end part of this curve at the parameter t.
//The new curve range is [start_of_original, t]
//t must be inside this parameter rage, else does nothing.
void MGCurve::trim_end(double t){
	MGInterval range = param_range();
	if(range<<t)
		return;
	range.set_high_point(t);
	limit(range);
}

//Trim the start part of this curve at the parameter t.
//The new curve range is [t,end_of_original]
//t must be inside this parameter rage, else does nothing.
void MGCurve::trim_start(double t){
	MGInterval range = param_range();
	if(range<<t)
		return;
	range.set_low_point(t);
	limit(range);
}

// Creates a curve that has weight.
//Returned object is a newed object. Use must delete it.
MGRLBRep* MGCL::convert_to_rational(const MGCurve& curve){
	long tid = curve.identify_type();
	switch(tid){
	case MGRLBREP_TID:
		return dynamic_cast<MGRLBRep*>(curve.copy_as_nurbs());
	case MGELLIPSE_TID:
		return new MGRLBRep(dynamic_cast<const MGEllipse&>(curve));
	case MGLBREP_TID:
		return new MGRLBRep(dynamic_cast<const MGLBRep&>(curve));
	default:
		{
			std::auto_ptr<MGRLBRep> r(new MGRLBRep(MGLBRep(curve)));
			if(r->order() < 4){
				r->change_order(4);
			}
			return r.release();
		}
	}
}
