/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/CompositeCurve.h"
#include "mg/Position_list.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/TrimmedCurve.h"
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

//
//Implements MGCompositeCurve Class.
//
//Defines MGCompositeCurve Class.
//MGCompositeCurve is a composite of other leaf curves.
//Assumedly they are connected as C0 continuity. However, MGCompositeCurve
//does not check their continuity, but only put first or last as the user says
// (in connect_to_end or in connect_to_start).
//Parameter ranges of the member curves are always continuous. param_s() of the
//1st curve to param_e() of the last form MGCompositeCurve's paramter range.
//number_of_curves() indicates the number of leaf curves(in other words,
//element curves to construct the MGCompositeCurve).

////////////Constructor/////////////

//Constructor of one curve.
//crv is a newed object pointer and MGCompositeCurve takes the ownership.
MGCompositeCurve::MGCompositeCurve(MGCurve* crv)
:MGCurve(*crv), m_composite(1,crv){
}

//Copy constructor.
MGCompositeCurve::MGCompositeCurve(const MGCompositeCurve& original)
:MGCurve(original),m_composite(original.m_composite.size()){
	const_iterator i=original.begin(), ie=original.end();
	for(size_t j=0; i!=ie; i++, j++)
		m_composite[j]=(**i).clone();
}

////////// Destructor //////////
MGCompositeCurve::~MGCompositeCurve(){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) delete *i;
}

////////// Operator overload(演算子多重定義) //////////

//Assignment.

//Assignment.
//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGCompositeCurve& MGCompositeCurve::operator=(const MGCompositeCurve& original){
	if(this==&original)
		return *this;

	MGCurve::operator=(original);
	iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		delete *i;
	size_t n=original.m_composite.size();
	m_composite.resize(n);
	const_iterator j=original.begin(), je=original.end();
	for(size_t k=0; j!=je; k++, j++)
		m_composite[k]=(**j).clone();
	return *this;
}
MGCompositeCurve& MGCompositeCurve::operator=(const MGGel& gel2){
	const MGCompositeCurve* gel2_is_this=dynamic_cast<const MGCompositeCurve*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線を平行移動して自身とする。
MGCompositeCurve& MGCompositeCurve::operator+= (const MGVector& v){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i)+=v;
	if(m_box) *m_box+=v;
	return *this;
}

//Update the curve by translation.
// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGCompositeCurve& MGCompositeCurve::operator-= (const MGVector& v){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i)-=v;
	if(m_box) *m_box-=v;
	return *this;
}

//Update the curve by multiplying scale.
// 与えられたスケールを曲線にかける。
MGCompositeCurve& MGCompositeCurve::operator*= (double scale){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i)*=scale;
	update_mark();
	return *this;
}

//Update the curve by transformation of matrix.
// 与えられた変換で直線の変換を行い自身の直線とする。
MGCompositeCurve& MGCompositeCurve::operator*= (const MGMatrix& mat){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i)*=mat;
	update_mark();
	return *this;
}

//Update the curve by transformation of transf.
// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGCompositeCurve& MGCompositeCurve::operator*= (const MGTransf& tr){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i)*=tr;
	update_mark();
	return *this;
}

//Transformation object construction
MGCompositeCurve MGCompositeCurve::operator+ (const MGVector& v)const{
	MGCompositeCurve cc(*this);
	cc+=v;
	return cc;
}
MGCompositeCurve operator+ (const MGVector& v, const MGCompositeCurve& cc){
	return cc+v;
}
MGCompositeCurve MGCompositeCurve::operator- (const MGVector& v) const{
	MGCompositeCurve cc(*this);
	cc-=v;
	return cc;
}
MGCompositeCurve MGCompositeCurve::operator* (double scale) const{
	MGCompositeCurve cc(*this);
	cc*=scale;
	return cc;
}
MGCompositeCurve operator* (double scale, const MGCompositeCurve& cc){
	return cc*scale;
}
MGCompositeCurve MGCompositeCurve::operator* (const MGMatrix& mat) const{
	MGCompositeCurve cc(*this);
	cc*=mat;
	return cc;
}
MGCompositeCurve MGCompositeCurve::operator* (const MGTransf& tr) const{
	MGCompositeCurve cc(*this);
	cc*=tr;
	return cc;
}

//Logical operator overload(論理演算子多重定義)
//Test if two curves are equal.
// 与曲線と自身が等しいかの比較判定を行う。

bool MGCompositeCurve::operator==(const MGTrimmedCurve& crv)const{
	return crv.operator==(*this);
}
bool MGCompositeCurve::operator==(const MGCompositeCurve& crv)const{
	size_t n=crv.m_composite.size();
	if(m_composite.size()!=n)
		return false;

	const_iterator i=begin(), ie=end(), j=crv.begin();
	for(; i!=ie; i++, j++)
		if((**i)!=(**j))
			return false;
	return true;
}
bool MGCompositeCurve::operator<(const MGCompositeCurve& gel2)const{
	size_t n1=number_of_curves(), n2=gel2.number_of_curves();
	if(n1==n2){
		if(n1==0)
			return false;
		return curve(0)<gel2.curve(0);
	}else
		return n1<n2;
}
bool MGCompositeCurve::operator==(const MGGel& gel2)const{
	const MGCompositeCurve* gel2_is_this=dynamic_cast<const MGCompositeCurve*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	else{
		const MGTrimmedCurve* gel2_is_compo=dynamic_cast<const MGTrimmedCurve*>(&gel2);
		if(gel2_is_compo)
			return operator==(*gel2_is_compo);

	}
	return false;
}
bool MGCompositeCurve::operator<(const MGGel& gel2)const{
	const MGCompositeCurve* gel2_is_this=dynamic_cast<const MGCompositeCurve*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

////////// Member Function ///////////

bool MGCompositeCurve::is_same_curve(const MGCurve& curve2)const{
	if(number_of_curves()!=1)
		return false;
	return m_composite[0]->operator==(curve2);
}

//Returns B-Rep Dimension.
//bdim of MGCompositeCurve is the sum of each member curves.
size_t MGCompositeCurve::bdim() const{
	const_iterator i=begin(), ie=end();
	size_t n=0;
	for(; i!=ie; i++) n+=(**i).bdim();
	return n;
}

//Return minimum box that includes the curve of parameter interval.
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
MGBox MGCompositeCurve::box_limitted(
	const MGInterval& rng	// Parameter Range of the curve.
	) const{
	size_t n=m_composite.size();
	if(!n) return MGBox();

	size_t i0,i1;
	if(rng.low().minus_infinite()) i0=0;
	else i0=find(rng.low_point());
	if(rng.high().plus_infinite()) i1=n-1;
	else i1=find(rng.high_point());
	if(i0==i1) return m_composite[i0]->box_limitted(rng);
	else{
		MGBox bx=m_composite[i0]->box_limitted(rng);
		bx|=m_composite[i1]->box_limitted(rng);
		for(size_t i=i0+1; i<i1; i++) bx|=m_composite[i]->box();
		return bx;	
	}
}

//Changing this object's space dimension.
MGCompositeCurve& MGCompositeCurve::change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i).change_dimension(sdim,start1,start2);
	update_mark();
	return *this;
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
void MGCompositeCurve::change_range(
	double t0,	//Parameter value for the start of original. 
	double t1)	//Parameter value for the end of original.
{
	size_t n=m_composite.size();
	if(!n) return;

	double oldsp=param_e()-param_s();
	double newsp=(t1-t0);
	double ratio=newsp/oldsp;
	double ts, te;
	size_t nm1=n-1;
	if(t0<t1){	//Case of no direction change.
		ts=t0;
		for(size_t i=0; i<nm1; i++){
			te=ts+m_composite[i]->param_span()*ratio;
			m_composite[i]->change_range(ts,te);
			ts=te;
		}
		m_composite[nm1]->change_range(ts,t1);
	}else{		//Case of direction change.
		ts=t1;
		size_t nhalf=n/2, i;
		for(i=0; i<nhalf; i++){//Change the curve ordering.
			MGCurve* crv=m_composite[i];
			size_t nm1mi=nm1-i;
			m_composite[i]=m_composite[nm1mi];
			m_composite[nm1mi]=crv;
		}
		for(i=0; i<nm1; i++){//Change each parameter range(and direction).
			double te=ts-ratio*(m_composite[i]->param_span());
			m_composite[i]->change_range(te,ts);
			ts=te;
		}
		m_composite[nm1]->change_range(t0,ts);
	}
	update_mark();
}

//Compute the closest point parameter value of this curve from a point.
double MGCompositeCurve::closest(const MGPosition& point) const{
	double param, dist=-1.;
	MGCParam_list list;
	const_iterator is=begin(), ie=end();
	for(const_iterator i=is; i!=ie; i++){
		double parami=(**i).closest(point);
		MGVector pdif=eval(parami)-point;
		if(i==is){
			param=parami;
			dist=pdif%pdif;
		}else{
			double disti=pdif%pdif;
			if(disti<dist){
				dist=disti;
				param=parami;
			}
		}
	}
	return param;
}

//Return minimum box that includes whole of the curve.
//曲線部分を囲むボックスを返す。
//Or'ed box of all the member curves.
MGBox* MGCompositeCurve::compute_box() const{
	MGBox* bx=new MGBox();
	if(!m_composite.size()) return bx;
	const_iterator i=begin(),	ie=end();
	for(; i!=ie; i++) (*bx) |=(*i)->box();
	return bx;	
}

//Connect the input curve to the end(to_end) or start(to_start) of this curve.
//(1) End(start) point of this curve is assumedly the same as the start(end) point
//of add_curve. However, connect_to_xxx does not check the continuity.
//(2) add_curve must be a newed object pointer and connect_to_end takes the ownership.
//(3) add_curve can be a MGCompositeCurve. In this case each curves in add_curve
//become members of this MGCompositeCurve.
//(4) When add_curve is a SurfCurve, or a TrimmedCurve, it is changed to non SurfCurve
//or non TrimmedCurve. Thus the original surface or curve can be deleted or modified.
//(5) When MGCompositeCurve was not an empty curve, the original part of
//the MGCompositeCurve's parameter range will not be changed by this connect_to_end.
//Instead add_curve's parameter range will be so modified that the magnitude of
//1st derivatives of the two curves are the same at the connecting point and
//the modified add_curve's parameter range is continuous to the original.
//(6) connect_to_end(start) will change add_curve's direction if necessary.
//Function's return value is the new parameter range of add_curve part after added.
//connect() connects to either start or end of this depending the distance
//of the start or end points of this curve to the start or end points of add_curve.
//connect_to_end() connects either start or end points of add_curve to the end to this curve.
//connect_to_start() connects either start or end points of add_curve to the start to this curve.
MGInterval MGCompositeCurve::connect(MGCurve* add_curve){
	size_t n=m_composite.size();
	if(n==0)
		return connect_to_end(add_curve);

	MGPosition P00=start_point(), P01=end_point();
	MGPosition P10=add_curve->start_point(), P11=add_curve->end_point();
	MGVector V0010=P00-P10, V0011=P00-P11;
	MGVector V0110=P01-P10, V0111=P01-P11;
	double dist0010=V0010%V0010, dist0011=V0011%V0011;
	double dist0110=V0110%V0110, dist0111=V0111%V0111;
	double error2=MGTolerance::wc_zero_sqr()*2.;
	if(dist0110<=error2 || dist0111<=error2)
		return connect_to_end(add_curve);
	if(n==1){
		if(dist0010<=error2 || dist0011<=error2){
			negate();
			return connect_to_end(add_curve);
		}
	}

	if((dist0010<=dist0110 && dist0010<=dist0111) ||
		(dist0011<=dist0110 && dist0011<=dist0111))
		return connect_to_start(add_curve);
	else
		return connect_to_end(add_curve);
}
MGInterval MGCompositeCurve::connect_to_end(MGCurve* add_curve){
	MGCompositeCurve* ccrv=dynamic_cast<MGCompositeCurve*>(add_curve);
	size_t n=m_composite.size();
	if(n==0){
		MGInterval rng=add_curve->param_range();
		if(ccrv){
			size_t n=ccrv->number_of_curves();
			for(size_t i=0; i<n; i++)
				connect_to_end(ccrv->release_front());
			delete add_curve;
		}else
			m_composite.push_back(add_curve);
		update_mark();
		return rng;
	}

	double tendold=param_e();
	double ts2=add_curve->param_s(), te2=add_curve->param_e();
	MGPosition Pend=eval(tendold);
	double a=eval(tendold,1).len();

	MGPosition Qstart=add_curve->eval(ts2);
	double b=add_curve->eval(ts2,1).len();
	MGPosition Qend=add_curve->eval(te2);
	double span2=te2-ts2;
	double tendnew=tendold;
	if(a<=MGTolerance::mach_zero())
		tendnew+=span2;
	else
		tendnew+=b*span2/a;

	bool reverse_direction=(Pend.distance(Qstart)<=Pend.distance(Qend));

	MGTrimmedCurve* tcrv=dynamic_cast<MGTrimmedCurve*>(add_curve);
	if(tcrv){//When trimmed curve.
		MGCurve* add_curve2=tcrv->copy_limitted();
		delete add_curve;
		add_curve=add_curve2;
	}
	if(ccrv){//When  add_curve is MGCompositeCurve.
		if(reverse_direction)
			ccrv->change_range(tendold, tendnew);
		else
			ccrv->change_range(tendnew, tendold);
		size_t n2=ccrv->m_composite.size();
		for(size_t i=0; i<n2; i++) m_composite.push_back(ccrv->m_composite[i]);
		ccrv->m_composite.clear();
		delete add_curve;
		update_mark();
		return MGInterval(tendold, tendnew);
	}
	MGSurfCurve* scrv=dynamic_cast<MGSurfCurve*>(add_curve);
	if(scrv){//Change to LBRep.
		MGLBRep* lb=new MGLBRep(*scrv);
		delete add_curve;
		add_curve=lb;
	}

	if(reverse_direction)
		add_curve->change_range(tendold, tendnew);
	else
		add_curve->change_range(tendnew, tendold);

	m_composite.push_back(add_curve);
	update_mark();
	return MGInterval(tendold, tendnew);
}
MGInterval MGCompositeCurve::connect_to_start(MGCurve* add_curve){
	MGCompositeCurve* ccrv=dynamic_cast<MGCompositeCurve*>(add_curve);
	size_t n=m_composite.size();
	if(n==0){
		MGInterval rng=add_curve->param_range();
		if(ccrv){
			size_t n=ccrv->number_of_curves();
			for(size_t i=0; i<n; i++)
				connect_to_end(ccrv->release_front());
			delete add_curve;
		}else
			m_composite.push_back(add_curve);
		update_mark();
		return rng;
	}
	double tstartold=param_s();
	double ts2=add_curve->param_s(), te2=add_curve->param_e();
	MGPosition Pstart=eval(tstartold);
	double a=eval(tstartold,1).len();

	MGPosition Qstart=add_curve->eval(ts2);
	MGPosition Qend=add_curve->eval(te2);
	double b=add_curve->eval(te2,1).len();
	double span2=te2-ts2;
	double tstartnew=tstartold;
	if(a<=MGTolerance::mach_zero()) tstartnew-=span2;
	else tstartnew-=b*span2/a;

	bool reverse_direction=(Pstart.distance(Qstart)<=Pstart.distance(Qend));

	MGTrimmedCurve* tcrv=dynamic_cast<MGTrimmedCurve*>(add_curve);
	if(tcrv){//When trimmed curve.
		MGCurve* add_curve2=tcrv->copy_limitted();
		delete add_curve;
		add_curve=add_curve2;
	}
	if(ccrv){//When  add_curve is MGCompositeCurve.
		if(reverse_direction)
			ccrv->change_range(tstartold,tstartnew);
		else
			ccrv->change_range(tstartnew,tstartold);
		size_t n2=ccrv->m_composite.size();
		for(int i=n2-1; i>=0; i--) m_composite.push_front(ccrv->m_composite[i]);
		ccrv->m_composite.clear();
		delete add_curve;
		update_mark();
		return MGInterval(tstartnew, tstartold);
	}
	MGSurfCurve* scrv=dynamic_cast<MGSurfCurve*>(add_curve);
	if(scrv){//Change to LBRep.
		MGLBRep* lb=new MGLBRep(*scrv);
		delete add_curve;
		add_curve=lb;
	}

	if(reverse_direction)
		add_curve->change_range(tstartold,tstartnew);
	else
		add_curve->change_range(tstartnew,tstartold);
	//cout<<"in connect_to_start::"<<(add_curve->eval(tstartold,1).len())<<endl;////////

	m_composite.push_front(add_curve);
	update_mark();
	return MGInterval(tstartnew, tstartold);
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGCompositeCurve& MGCompositeCurve::coordinate_exchange(size_t i, size_t j){
	iterator k=begin(), ke=end();
	for(; k!=ke; k++) (**k).coordinate_exchange(i,j);
	update_mark();
	return *this;
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGCompositeCurve* MGCompositeCurve::clone() const{
	return new MGCompositeCurve(*this);
}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGCompositeCurve::copy_as_nurbs() const{
	size_t n=number_of_curves();
	std::vector<MGCurve*> crvs(n);
	bool homo=true;//true means all of the constituent curves are MGLBRep, or MGRLBRep.
	crvs[0]=curve(0).copy_as_nurbs();
	MGLBRep* lb0=dynamic_cast<MGLBRep*>(crvs[0]);
	MGRLBRep* rlb0=dynamic_cast<MGRLBRep*>(crvs[0]);

	//Store newed NURBS in crvs.
	size_t j;
	for(j=1; j<n; j++){
		crvs[j]=curve(j).copy_as_nurbs();
		if(homo){
			if(dynamic_cast<MGLBRep*>(crvs[j])){
				if(rlb0)
					homo=false;
			}else{
				if(lb0)
					homo=false;
			}
		}
	}

	MGCurve* new_curve;
	//int which; double ratio;
	if(homo){
		if(lb0){//If all of the curves are MGLBRep.
			MGLBRep *lbj;
			for(j=1; j<n; j++){
				lbj=dynamic_cast<MGLBRep*>(crvs[j]);
				//int continuity=lb0->continuity(*lbj, which,ratio);
				lb0->connect(0,2,*lbj);
				delete crvs[j];
			}
			new_curve=lb0;
		}else{//If all of the curves are MGRLBRep.
			for(j=1; j<n; j++){
				MGRLBRep* rlbj=dynamic_cast<MGRLBRep*>(crvs[j]);
				//int continuity=rlb0->continuity(*rlbj, which,ratio);
				rlb0->connect(0,2,*rlbj);
				delete crvs[j];
			}
			new_curve=rlb0;
		}
	}else{//MGRLBRep and MGLBRep are mixed. MGLBRep will be changed to MGRLBRep.
		MGRLBRep* rlb;
		if(rlb0){
			rlb=rlb0;
		}else{
			rlb=new MGRLBRep(*lb0);
			delete crvs[0];
		}
		for(j=1; j<n; j++){
			MGRLBRep *rlbj=dynamic_cast<MGRLBRep*>(crvs[j]);
			if(!rlbj){
				rlbj=new MGRLBRep(*(dynamic_cast<MGLBRep*>(crvs[j])));
				delete crvs[j];
				crvs[j]=rlbj;
			}
			//int continuity=rlb->continuity(*rlbj, which,ratio);
			rlb->connect(0,2,*rlbj);
			delete crvs[j];
		}
		new_curve=rlb;
	}
	return new_curve;
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGCompositeCurve* MGCompositeCurve::copy_change_dimension(
	size_t sdim,	// new space dimension
	size_t start1, 	// Destination order of new line.
	size_t start2	// Source order of this line.
)const{ 		
	
	MGCompositeCurve* ccrv=new MGCompositeCurve;
	const_iterator k=begin(), ke=end();
	for(; k!=ke; k++)
		ccrv->m_composite.push_back(
		(**k).copy_change_dimension(sdim,start1,start2));
	return ccrv;
}

//Compute curvilinear integral of the 1st two coordinates from parameter t0
//to t1.
//This integral can be used to compute area sorrounded by the curve.
//The sum of all the member curve's curvilinear_integral.
double MGCompositeCurve::curvilinear_integral(double t0, double t1) const{
	if(m_composite.size()==0) return 0.;

	double s0,s1;
	if(t0>t1){ s0=t1; s1=t0;}
	else{ s0=t0; s1=t1;}
	size_t i0=find(s0), i1=find(s1);
	if(i0==i1) return m_composite[i0]->curvilinear_integral(s0,s1);
	else{
		double integral=m_composite[i0]->
			curvilinear_integral(s0, m_composite[i0]->param_e());
		for(size_t i=i0+1; i<i1; i++)
			integral+=m_composite[i]->curvilinear_integral();
		integral+=m_composite[i1]->
			curvilinear_integral(m_composite[i1]->param_s(), s1);
		if(t0>t1) integral*=-1.;
		return integral;
	}
}

//Compute curvilinear integral of the 1st two coordinates.
//(All the parameter range of the curve.)
//This integral can be used to compute area sorrounded by the curve.
//The sum of all the member curve's curvilinear_integral.
double MGCompositeCurve::curvilinear_integral() const{
	double integral=0.;
	const_iterator i=begin(),	ie=end();
	for(; i!=ie; i++)
		integral+=(**i).curvilinear_integral();
	return integral;
}

void MGCompositeCurve::display_break_points()const{
	const_iterator i=begin(),	ie=end();
	for(; i!=ie; i++)
		(**i).display_break_points();
}
void MGCompositeCurve::display_control_polygon()const{
	const_iterator i=begin(),	ie=end();
	for(; i!=ie; i++)
		(**i).display_control_polygon();
}
void MGCompositeCurve::display_curvatures(
	double	scale,	//scaling of the graph.
	int		density,//densitiy of the graph.
	bool	use_radius//true:radius display, false:curvature display.
)const{
	const_iterator i=begin(),	ie=end();
	for(; i!=ie; i++)
		(**i).display_curvatures(scale,density,use_radius);
}

//Divide this curve at the designated knot multiplicity point.
//Function's return value is the number of the curves after divided.
int MGCompositeCurve::divide_multi(
	MGPvector<MGCurve>& crv_list,	//divided curves will be appended.
	int multiplicity	//designates the multiplicity of the knot to divide at.
						//When multiplicity<=0, order()-1 is assumed.
						//When multiplicity>=order(), order() is assumed.
) const{
	const_iterator i=begin(),	ie=end();
	int num=0;
	for(; i!=ie; i++)
		num+=(**i).divide_multi(crv_list,multiplicity);
	return num;
}

//Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGCompositeCurve::eval(
	double t,			// Parameter value.
	size_t nderiv,		// Order of Derivative.
	int left			//Left continuous(left=true)
						//or right continuous(left=false).
)const{
	size_t ncrv=m_composite.size();
	if(!ncrv)
		return MGVector();

	double pspan=param_span();
	size_t i=find(t);
	const MGCurve& curvei=*(m_composite[i]);
	if(left){
		if(0<i && MGREqual_base(t,curvei.param_s(),pspan))
			return m_composite[i-1]->eval(t,nderiv,left);
	}else{
		if(i<ncrv-1 && MGREqual_base(t,curvei.param_e(),pspan))
			return m_composite[i+1]->eval(t,nderiv,left);
	}
	return curvei.eval(t,nderiv,left);
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGCompositeCurve::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	size_t n=number_of_curves();
	if(!n)
		return;
	size_t i=start ? 0:n-1;
	curve(i).extend(length,start);
}

//Find which curve parameter range input t belongs to.
//Function's return value m is the id of m_composite which satisfy
//m_composite[m]<=t<m_composite[m+1];
//m_composite.size() should be >=1, or find will abort.
size_t MGCompositeCurve::find(double t, bool right_continuous) const{
	size_t high=m_composite.size(), low=0;
	if(high<=1)
		return 0;

	assert(m_composite.front()->param_s()<m_composite.back()->param_e());
	size_t m=0;
	if(t<=m_composite.front()->param_e())
		m=low;
	else if(t>=m_composite.back()->param_s())
		m=high-1;
	else{
		//Perform the binary search.
		while((high-low)>1){
			m=(low+high)/2;
			const MGCurve& crv=*(m_composite[m]);
			if(t<crv.param_s()){
				high=m; continue;
			}
			if(crv.param_e()<=t){
				low=m; continue;
			}
			break;
		}
	}
	double perror=param_error();
	const MGCurve& crv=*(m_composite[m]);
	if(right_continuous && m<high-1){
		if(t>crv.param_e()-perror)
			m++;
	}else if((!right_continuous) && m>0){
		if(t<crv.param_s()+perror)
			m--;
	}
	return m;
}

//Provide divide number of curve span for function intersect.
//This will never be used since each member's function will be used.
size_t MGCompositeCurve::intersect_dnum() const{
	size_t num=0;
	const_iterator i=begin(),	ie=end();
	for(; i!=ie; i++)
		num+= (**i).intersect_dnum();
	return num;
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect_1by1(const MGCurve& crv2) const{
	MGCCisect_list list(this, &crv2);
	MGBox bx=box(); bx&=crv2.box();
	if(bx.empty())
		return list;

	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		list.append((**i).isect(crv2));
	return list;
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGCurve& crv2) const{
	return isect_1by1(crv2);
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGStraight& crv2) const{
	return isect_1by1(crv2);
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list MGCompositeCurve::isect(const MGRLBRep& curve2)const{
	return isect_1by1(curve2);
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGEllipse& crv2) const{
	return isect_1by1(crv2);
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list MGCompositeCurve::isect(const MGLBRep& curve2)const{
	return isect_1by1(curve2);
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGSurfCurve& crv2) const{
	return isect_1by1(crv2);
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGBSumCurve& crv2) const{
	return isect_1by1(crv2);
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGTrimmedCurve& crv2) const{
	return isect_1by1(crv2);
}

//Intersection of CompositeCurve and another curve crv2.
MGCCisect_list MGCompositeCurve::isect(const MGCompositeCurve& crv2) const{
	return isect_1by1(crv2);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGSurface& surf) const{
	return isect_1by1(surf);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGPlane& surf) const{
	return isect_1by1(surf);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGSphere& surf) const{
	return isect_1by1(surf);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGCylinder& surf) const{
	return isect_1by1(surf);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGSBRep& surf) const{
	return isect_1by1(surf);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGRSBRep& surf) const{
	return isect_1by1(surf);
}

//Intersection with a Surface
MGCSisect_list MGCompositeCurve::isect(const MGBSumSurf& surf) const{
	return isect_1by1(surf);
}

//Intersection with a surface. Obtains the isects by computing isects
//of each elemet curve.
MGCSisect_list MGCompositeCurve::isect_1by1(const MGSurface& surf)const{
	MGCSisect_list list(this, &surf);
	MGBox bx=box(); bx&=surf.box();
	if(bx.empty())
		return list;

	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		list.append((**i).isect(surf));
	return list;
}

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
MGCParam_list MGCompositeCurve::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
) const{
	MGCParam_list list(this);
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		list.append((**i).isect_1D(f,coordinate));
	return list;
}

//Test if this cure is planar or not.
//MGPlane expression will be output to plane if this is planar.
//Function's return value is true if planar.
bool MGCompositeCurve::is_planar(MGPlane& plane)const{
	const_iterator i=begin(), ie=end();
	MGStraight line;
	MGPosition point;
	bool plane_is_obtained=false;
	for(; i!=ie; i++){
		if(!(**i).is_planar(plane))
			return false;

		const MGStraight* sl=dynamic_cast<const MGStraight*>(*i);
		if(sl) continue;
		const MGLBRep* lb=dynamic_cast<const MGLBRep*>(*i);
		if(lb){
			int plkind=lb->planar(plane,line,point);
			if(plkind==3){
				plane_is_obtained=true;
				break;
			}else
				continue;
		}
		const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(*i);
		if(rlb){
			int plkind=rlb->planar(plane,line,point);
			if(plkind==3){
				plane_is_obtained=true;
				break;
			}else
				continue;
		}
	}
	if(number_of_curves()<=1)
		return true;

	i=begin();
	if(!plane_is_obtained){
		const MGCurve& c0=*(*i++);
		const MGCurve& c1=*(*i++);
		MGVector udir=c0.eval_deriv(c0.param_e());
		MGVector vdir=c1.eval_deriv(c1.param_s());
		plane=MGPlane(udir,vdir,c0.eval(c0.param_e()));
	}
	MGPlane plane2;
	MGPosition uv;
	for(; i!=ie; i++){
		const MGStraight* sl=dynamic_cast<const MGStraight*>(*i);
		if(sl){
			if(!plane.on(*sl))
				return false;
			continue;
		}
		const MGLBRep* lb=dynamic_cast<const MGLBRep*>(*i);
		if(lb){
			int plkind=lb->planar(plane2,line,point);
			if(plkind==1){
				if(!plane.on(point,uv))
					return false;
				continue;
			}else if(plkind==2){
				if(!plane.on(line))
					return false;
				continue;
			}else{
				if(plane2!=plane)
					return false;
				continue;
			}
		}
		const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(*i);
		if(rlb){
			int plkind=rlb->planar(plane2,line,point);
			if(plkind==1){
				if(!plane.on(point,uv))
					return false;
				continue;
			}else if(plkind==2){
				if(!plane.on(line))
					return false;
				continue;
			}else{
				if(plane2!=plane)
					return false;
				continue;
			}
		}
	}
	return true;
}

//Access to i-th element of knot.
double MGCompositeCurve::knot(size_t i) const{
	if(i<order()) return param_s();
	if(i>=bdim()) return param_e();
	if(number_of_curves()==1) return m_composite[0]->knot(i);
//	assert(false);
	return param_s();
}

//Returns the knot vector of the curve.
//This should not be used.
const MGKnotVector& MGCompositeCurve::knot_vector() const{
	if(number_of_curves()==1)
		return m_composite[0]->knot_vector();
//	assert(false);
	return mgNULL_KNOT_VECTOR;
}

//Cmpute curve length of the interval.
//If t1 is greater than t2, return negative value.
// パラメータが昇順で与えられたときは正値、降順のときは負値を返す。
//The sum of all the member curve's length.
double MGCompositeCurve::length(double t0, double t1) const{
	if(m_composite.size()==0) return 0.;

	double s0, s1;
	if(t0<t1){s0=t0; s1=t1;} else{ s0=t1; s1=t0;}
	size_t i0=find(s0), i1=find(s1);
	if(i0==i1) return m_composite[i0]->length(t0,t1);
	else{
		double len=m_composite[i0]->length(s0, m_composite[i0]->param_e());
		for(size_t i=i0+1; i<i1; i++)
			len+=m_composite[i]->length();
		len+=m_composite[i1]->length(m_composite[i1]->param_s(), s1);
		if(t1<t0) len*=-1.;
		return len;
	}
}

//Update this by limiting the parameter range of the curve.
// 自身に指定したパラメータ範囲のｌｉｍｉｔをつける。
MGCompositeCurve& MGCompositeCurve::limit(const MGInterval& rng){
	if(m_composite.size()==0) return *this;

	double t0=rng.low_point(), t1=rng.high_point();
	size_t i0,i1;
	if(rng>param_s()){
		i0=find(t0);
		for(size_t i=0; i<i0; i++){
			MGCurve* c1st=m_composite.front();
			m_composite.pop_front();
			delete c1st;
		}
	}
	if(rng<param_e()){
		size_t n=m_composite.size();
		i1=find(t1,false);
		for(size_t i=i1+1; i<n; i++){
			MGCurve* clast=m_composite.back();
			m_composite.pop_back();
			delete clast;
		}
	}
	m_composite.front()->limit(rng);
	if(m_composite.size()>1)
		m_composite.back()->limit(rng);
	update_mark();
	return *this;
}

//Negate the curve direction(曲線の方向を反転する)
void MGCompositeCurve::negate(){
	size_t n=m_composite.size();
	if(!n) return;

	size_t i, nm1=n-1;
	double ts=param_s(), teOrigin=param_e();

	size_t nhalf=n/2;
	for(i=0; i<nhalf; i++){//Reverse the curve ordering.
		MGCurve* crv=m_composite[i];
		size_t nm1mi=nm1-i;
		m_composite[i]=m_composite[nm1mi];
		m_composite[nm1mi]=crv;
	}

	for(i=0; i<nm1; i++){//Change each parameter range(and direction).
		double te=ts+m_composite[i]->param_span();
		m_composite[i]->change_range(te,ts);
		ts=te;
	}
	m_composite[nm1]->change_range(teOrigin,ts);
}

//Obtain parameter value if this curve is negated by "negate()".
double MGCompositeCurve::negate_param(double t)const{
	int n=m_composite.size();
	if(!n) return t;

	int i=find(t);
	double ts=param_s();
	for(int j=n-1; j>i; j--)
		ts+=(m_composite[j]->param_e() - m_composite[j]->param_s());
	return ts+(m_composite[i]->param_e()-t);
}

//Test if given point is on the curve or not. If yes, return parameter
//value of the curve. Even if not, return nearest point's parameter.
// 指定点が自身上にあるかを調べる。曲線上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the curve,
// and 0 if the point is not on the curve.
bool MGCompositeCurve::on(
	const MGPosition& P,	//Point(指定点)
	double& t				//Parameter of the curve(パラメータ)
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((**i).on(P,t))
			return true;
	}
	return false;
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGCompositeCurve::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
)const{
	const_iterator i=begin(), ie=end();
	MGCompositeCurve* compo=new MGCompositeCurve();
	for(; i!=ie; i++){
		std::auto_ptr<MGCurve> crv=(*i)->oneD(g);
		compo->m_composite.push_back(crv.release());
	}
	return std::auto_ptr<MGCurve>(compo);
}

//Returns the order.
//Returns the maximum order among the curves.
unsigned MGCompositeCurve::order() const{
	const_iterator i=begin(), ie=end();
	size_t n=0;
	for(; i!=ie; i++){
		size_t n1=(**i).order();
		if(n<n1) n=n1;
	}
	return n;
}

//Return ending parameter value.
double MGCompositeCurve::param_e() const{
	if(m_composite.size()==0) return 0.;
	return m_composite.back()->param_e();
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance.
double MGCompositeCurve::param_normalize(double t) const{
	if(m_composite.size()==0) return 0.;

	size_t i=find(t);
	return m_composite[i]->param_normalize(t);	
}

// Return starting parameter value.
double MGCompositeCurve::param_s() const{
	if(m_composite.size()==0) return 0.;
	return m_composite.front()->param_s();
}
	
//Compute part of this curve from parameter t0 to t1.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGCurve* MGCompositeCurve::part(
	double t0, double t1,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
) const{
	assert(t0<=t1);
	if(t0<=param_s() && param_e()<=t1)
		return clone();

	size_t n=m_composite.size();
	if(n==0)
		return clone();

	size_t i0=find(t0), i1=find(t1,false);
	if((i1-i0)<=0)
		return m_composite[i0]->part(t0,t1,multiple);

	MGCompositeCurve* ccrv=new MGCompositeCurve();
	ccrv->m_composite.push_back(m_composite[i0]->part(t0,t1,multiple));
	for(size_t i=i0+1; i<i1; i++){
		ccrv->m_composite.push_back(m_composite[i]->clone());
	}
	ccrv->m_composite.push_back(m_composite[i1]->part(t0,t1,multiple));
	return ccrv;
}

//Compute all foot points of the perpendicular line from point to
//the curve.
// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
// すべて求める。
MGCParam_list MGCompositeCurve::perps(
	const MGPosition& P		//Point(指定点)
) const{
	MGCParam_list list;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		list.append((**i).perps(P));
	return list;
}

MGPosition_list MGCompositeCurve::perps(const MGSurfCurve& curve2)const{
	return perps_1by1(curve2);	
}

//Return perpendicular point from a point P,
//given guess starting paramter values.
//Function's return value is:
//   perp_guess=true if perpendicular points obtained,
//   perp_guess=false if perpendicular points not obtained,
int MGCompositeCurve::perp_guess(
	double t0, double t1,	//parameter range of this.
							//(t0>=t1) indicates no range specified.
	const MGPosition& P,	//Point(指定点)
	double tg,				//Guess parameter values of this curve.
	double& t				//Output parameter
)const{
	size_t cnum=find(tg);
	return m_composite[cnum]->perp_guess(t0,t1,P,tg,t);
}

//Return perpendicular points of two curves,
//given guess starting paramter values.
//Function's return value is:
//   perp_guess=true if perpendicular points obtained,
//   perp_guess=false if perpendicular points not obtained,
int MGCompositeCurve::perp_guess(
	double s0, double s1,		//parameter range of this.
					//When s0>=s1, no limit for this parameter range.
	const MGCurve& curve2,		//2nd curve.
	double t0, double t1,		//parameter range of curve2.
					//When t0>=t1, no limit for curve2 parameter range.
	double sg, double tg,		//Guess parameter values of the two curves
	//sg: this curve's parameter, tg:curve2's parameter.
	MGPosition& st				//perpendicular points' parameter values
								//will be output.
	//st(0): this curve's parameter, st(1):curve2's parameter.
)const{
	size_t cnum=find(sg);
	return m_composite[cnum]->perp_guess(s0,s1,curve2,t0,t1,sg,tg,st);
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list MGCompositeCurve::perps_1by1(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		list.append(*this,crv2,(**i).perps(crv2));
	return list;
}

//Release the pointer of the last curve.
//Returned will be the released MGCurve pointer.
MGCurve* MGCompositeCurve::release_back(){
	iterator i=end(); i--;
	MGCurve* crv=*i;
	m_composite.erase(i);
	return crv;
}

//Release the pointer of the 1st curve.
//Returned will be the released MGCurve pointer.
MGCurve* MGCompositeCurve::release_front(){
	iterator i=begin();
	MGCurve* crv=*i;
	m_composite.erase(i);
	return crv;
}

//ノット削除関数(B表現曲線のみ)
//トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
//Remove redundant knot, and reduce the b-rep dimension.
//The tolerance used is MGTolerance::line_zero().
void MGCompositeCurve::remove_knot(){
	iterator i=begin(), ie=end();
	for(; i!=ie; i++) (**i).remove_knot();
	update_mark();
}

//Return the space dimension. It is the maximum space dimension of the
//member curves.
size_t MGCompositeCurve::sdim() const{
	const_iterator i=begin(),ie=end();
	size_t n=0;
	for(; i!=ie; i++){
		size_t n1=(**i).sdim();
		if(n<n1) n=n1;
	}
	return n;
}

//Unlimit parameter range of the curve(limitをはずす)
MGCurve& MGCompositeCurve::unlimit(){
	size_t n=m_composite.size();
	if(!n) return *this;

	m_composite.front()->unlimit_start();
	m_composite.back()->unlimit_end();
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGCompositeCurve::unlimit_end(){
	size_t n=m_composite.size();
	if(!n) return *this;

	m_composite.back()->unlimit_end();
	update_mark();
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGCompositeCurve::unlimit_start(){
	size_t n=m_composite.size();
	if(!n) return *this;

	m_composite.front()->unlimit_start();
	update_mark();
	return *this;
}
