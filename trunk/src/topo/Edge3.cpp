/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Bisection.h"
#include "topo/Edge.h"
#include "topo/Face.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

///////////////////////////////////////////////////
//MGEdge::compute_continuity//////////////////////

class MG2CrvBisect:public MGBisection{
public:

	MG2CrvBisect(
		const MGCurve& curve1,
		const MGCurve& curve2
	);

	//evaluate the target value to find the solution.
	virtual double eval(
		double t1,	//parameter value of m_curve1.
		double t2	//parameter value of m_curve2.
	)const=0;

	//This is used in compare_replace() when new_value is evaluated. necessary_to_replace()
	//decides if the new value is more suitable solution than the old one.
	//When this returns true, it means the new value is more suitable solution.
	virtual bool necessary_to_replace(double new_value)const{return new_value>m_value;};

	//Prepare for the function compare_replace() by setting the initial value.
	virtual void set_initial_t(double t);

	//compare with the previous function value(the initial value is set
	//by set_initial_t) and replace t with the previous one if necessary.
	//The function's return value is the new parameter value.
	virtual double compare_replace(
		double t,		//parameter value to compare at of m_curve1.
		bool& replaced	//true will be returned if the t is the new solution candidate value.
	);

	//get the result of the solution.
	void get_result(double& curve1_t, double& curve2_t, double& value)const{
		curve1_t=m_curve1_t; curve2_t=m_curve2_t; value=m_value;
	}

	const MGCurve& m_curve1;
	double m_curve1_t;//first curve's parameter of m_value.

	const MGCurve& m_curve2;
	double m_curve2_t;//second curve's parameter of m_value.

	double m_value;
};

MG2CrvBisect::MG2CrvBisect(
	const MGCurve& curve1,
	const MGCurve& curve2
):MGBisection(curve1.param_s(),curve1.param_e()),
m_value(-1.),m_curve1(curve1), m_curve2(curve2){;}

//Prepare for the function compare_replace() by setting the initial value.
void MG2CrvBisect::set_initial_t(double t){
	m_curve1_t=t;
	MGPosition P1=m_curve1.eval(m_curve1_t);
	m_curve2_t=m_curve2.closest(P1);
	m_value=eval(m_curve1_t, m_curve2_t);
}

//compare with the previous function value(the initial value is set
//by set_initial_t) and replace t with the previous one if necessary.
//The function's return value is the new parameter value.
double MG2CrvBisect::compare_replace(
	double t,		//parameter value to compare at of m_curve1.
	bool& replaced	//true will be returned if the t is the new solution candidate value.
){
	MGPosition P1=m_curve1.eval(t);
	double t2=m_curve2.closest(P1);
	double value=eval(t,t2);
	replaced=necessary_to_replace(value);
	if(replaced){
		m_value=value;
		m_curve1_t=t;
		m_curve2_t=t2;
	}
	return m_curve1_t;
}

class MGEdgeMaxDistance:public MG2CrvBisect{
public:

	MGEdgeMaxDistance(
		const MGCurve& curve1,
		const MGCurve& curve2
	);

	//evaluate the distance of curve1 and 2.
	double eval(
		double t1,	//parameter value of m_curve1.
		double t2	//parameter value of m_curve2.
	)const;
};

MGEdgeMaxDistance::MGEdgeMaxDistance(
	const MGCurve& curve1,
	const MGCurve& curve2
):MG2CrvBisect(curve1,curve2){;}

//evaluate the distance of curve1 and 2.
double MGEdgeMaxDistance::eval(
	double t1,	//parameter value of m_curve1.
	double t2	//parameter value of m_curve2.
)const{
	MGPosition P1=m_curve1.eval(t1);
	MGPosition P2=m_curve2.eval(t2);
	return P1.distance(P2);
}

class MGEdgeMinDistance:public MGEdgeMaxDistance{
public:

	MGEdgeMinDistance(
		const MGCurve& curve1,
		const MGCurve& curve2
	):MGEdgeMaxDistance(curve1,curve2){;};

	//This is used in compare_replace() when new_value is evaluated. necessary_to_replace()
	//decides if the new value is more suitable solution than the old one.
	//When this returns true, it means the new value is more suitable solution.
	virtual bool necessary_to_replace(double new_distance)const{return new_distance<m_value;};
};

class MGEdgeTanDifference:public MG2CrvBisect{
public:

	MGEdgeTanDifference(
		const MGCurve& curve1,
		const MGCurve& curve2
	);

	//evaluate the tangent difference of curve1 at t1 and curve2 at t2.
	double eval(
		double t1,	//parameter value of m_curve1.
		double t2	//parameter value of m_curve2.
	)const;
};

MGEdgeTanDifference::MGEdgeTanDifference(
	const MGCurve& curve1,
	const MGCurve& curve2
):MG2CrvBisect(curve1,curve2){;}

double MGEdgeTanDifference::eval(
	double t1,	//parameter value of m_curve1.
	double t2	//parameter value of m_curve2.
)const{
	MGVector Tan1=m_curve1.eval(t1,1);
	MGVector Tan2=m_curve2.eval(t2,1);
	if(Tan1%Tan2<0.)
		Tan2*=-1.;
	 return Tan1.angle(Tan2);
}

class MGEdgeNormalDifference:public MG2CrvBisect{
public:

	MGEdgeNormalDifference(
		const MGEdge& edge1,	//Parameter edge1.
		const MGCurve& curve1,	//edge1's binder edge's curve.

		const MGEdge& edge2,	//Parameter edge2.
		const MGCurve& curve2	//edge2's binder edge's curve.
	);

	//evaluate the normal difference of face1 at t1 and face2 at t2.
	double eval(
		double t1,	//parameter value of m_curve1.
		double t2	//parameter value of m_curve2.
	)const;

	const MGFace& m_face1;
	const MGEdge& m_edge1;	//Parameter edge1.

	const MGFace& m_face2;
	const MGEdge& m_edge2;	//Parameter edge2.
};

MGEdgeNormalDifference::MGEdgeNormalDifference(
	const MGEdge& edge1,	//Parameter edge1.
	const MGCurve& curve1,	//Binder edge1's curve.

	const MGEdge& edge2,	//Parameter edge2.
	const MGCurve& curve2	//Binder edge2's curve.
):MG2CrvBisect(curve1,curve2),
m_face1(*edge1.face()),m_edge1(edge1),
m_face2(*edge2.face()),m_edge2(edge2){;}

//evaluate the normal difference of face1 at t1 and face2 at t2.
double MGEdgeNormalDifference::eval(
	double t1,	//parameter value of m_curve1.
	double t2	//parameter value of m_curve2.
)const{
	//max ndif(normal difference) info.
	double t1p=m_edge1.param_pcell(t1);
	MGVector N1=m_face1.normal(m_edge1.eval(t1p));

	double t2p=m_edge2.param_pcell(t2);
	MGVector N2=m_face2.normal(m_edge2.eval(t2p));
	if(N1%N2<0.)
		N2*=-1.;
	return N1.angle(N2);
}

//Exclusive function for MGEdge::compute_continuity.
//Obtain the comparison span of the two curves, crv1 and crv2.
void get_comparison_span(
	const MGCurve& crv1,
	const MGCurve& crv2,
	MGInterval& sspan,//span of crv1 will be output.
	MGInterval& tspan //span of crv2 will be output.
){
//Here, sx means crv1's parameter value, and tx means crv2's.

	double sS=crv1.param_s(), sE=crv1.param_e();
	double tS=crv2.param_s(), tE=crv2.param_e();
	bool c1closed=crv1.is_closed();
	bool c2closed=crv2.is_closed();
	if(c1closed&&c2closed){
		sspan=crv1.param_range();
		tspan=crv2.param_range();
		return;
	}

	if(c1closed && !c2closed){
		sspan=MGInterval(sS);//null interval.
		tspan=crv2.param_range();
		return;
	}
	if(!c1closed && c2closed){
		sspan=crv1.param_range();
		tspan=MGInterval(tS);//null interval.
		return;
	}

	MGPosition P1S=crv1.start_point(), P1E=crv1.end_point();
	double t0=crv2.closest(P1S), t1=crv2.closest(P1E);
	double DsSt0=P1S.distance(crv2.eval(t0)), DsEt1=P1E.distance(crv2.eval(t1));
	//(sS, t0, DsSt0), (sE, t1, DsEt1) make pairs.

	MGPosition P2S=crv2.start_point(), P2E=crv2.end_point();
	double s0=crv1.closest(P2S), s1=crv1.closest(P2E);
	double Ds0tS=P2S.distance(crv1.eval(s0)), Ds1tE=P2E.distance(crv1.eval(s1));
	//(s0,tS,Ds0tS), (s1,tE,Ds1tE) make pairs.

	//Define the direction from the smallest length parameter position.
	bool same_direction=crv1.has_same_direction_at(s1,crv2,tE);//when Ds1tE is the smallest.
	if(DsSt0<=DsEt1){
		if(DsSt0<=Ds0tS){
			if(DsSt0<=Ds1tE)
				same_direction=crv1.has_same_direction_at(sS,crv2,t0);//DsSt0 is the smallest.
		}else{
			if(Ds0tS<=Ds1tE)
				same_direction=crv1.has_same_direction_at(s0,crv2,tS);//Ds0tS is the smallest.
		}
	}else{
		if(DsEt1<=Ds0tS){
			if(DsEt1<=Ds1tE)
				same_direction=crv1.has_same_direction_at(sE,crv2,t1);//DsEt1 is the smallest.
		}else{
			if(Ds0tS<=Ds1tE)
				same_direction=crv1.has_same_direction_at(s0,crv2,tS);//Ds0tS is the smallest.
		}
	}

	if(same_direction){
		if(DsSt0>Ds0tS){
			tspan=MGInterval(tS,t1);
			sspan=MGInterval(s0,(s1<s0 ? sE:s1));
		}else{
			sspan=MGInterval(sS,s1);
			tspan=MGInterval(t0,(t1<t0 ? tE:t1));
		}
	}else{
		if(DsSt0>Ds1tE){
			tspan=MGInterval(t1,tE);
			sspan=MGInterval(s1,(s0>s1 ? s0:sE));
		}else{
			sspan=MGInterval(sS,s0);
			tspan=MGInterval((t0>t1 ? t1:tS),t0);
		}
	}
}

//Compute the continuities between this edge(edge1) and the edge2.
//This edge and edge2 must be parameter edges of each face.
//In distance, tangent, and normal, the following output will be set:
//distance[0-6] as:
//	[0] edge1's binder edge's curve parameter that has the maximum distance with edge2.
//	[1] edge2's binder edge's curve parameter that has the maximum distance with edge1.
//  [2] the evaluated maximum distance between edge1 and edge2 at distance[0] and [1]
//	[3] edge1's binder edge's curve parameter that has the minimum distance with edge2.
//	[4] edge2's binder edge's curve parameter that has the minimum distance with edge1.
//  [5] the evaluated minimum distance between edge1 and edge2 at distance[3] and [4]
//	[6] mean distance between edge1 and edge2.
//tangent[0-3] as:
//	[0] edge1's binder edge's curve parameter that has the maximum tangent difference with edge2.
//	[1] edge2's binder edge's curve parameter that has the maximum tangent difference with edge1.
//  [2] the evaluated maximum tangent difference between edge1 and edge2 at tangent[0] and [1].
//	[3] mean tangent difference between edge1 and edge2.
//normal[0-3] as:
//	[0] edge1's binder edge's curve parameter that has the maximum normal difference with edge2.
//	[1] edge2's binder edge's curve parameter that has the maximum normal difference with edge1.
//  [2] the evaluated maximum normal difference between edge1 and edge2 at normal[0] and [1].
//	[3] mean normal difference between edge1 and edge2.
void MGEdge::compute_continuity(
	const MGEdge& edge2,	//Parameter edges of the 2nd face.
	double distance[7],
	double tangent[4],
	double normal[4]
)const{
	MGEdge* be1=make_binder_with_curve();
	std::auto_ptr<MGCurve> crv1(be1->curve_limitted());//world curve repr of edge1.

	MGEdge* be2=edge2.make_binder_with_curve();
	std::auto_ptr<MGCurve> crv2(be2->curve_limitted());//world curve repr of edge1.

	MGInterval sspan, tspan;
	get_comparison_span(*crv1,*crv2,sspan,tspan);

	int nspoint=crv1->offset_div_num(sspan);
	int ntpoint=crv2->offset_div_num(tspan);
	if(nspoint>ntpoint){
		compute_continuity2(sspan,nspoint,edge2,distance,tangent,normal);
	}else{
		edge2.compute_continuity2(tspan,ntpoint,*this,distance,tangent,normal);
		double tmp=distance[0];distance[0]=distance[1]; distance[1]=tmp;
		tmp=distance[3];distance[3]=distance[4]; distance[4]=tmp;
		tmp=tangent[0];tangent[0]=tangent[1]; tangent[1]=tmp;
		tmp=normal[0];normal[0]=normal[1]; normal[1]=tmp;
	}
}

//Compute continuity, given the evaluation interval and the division number.
//This is a parameter edge that has the star face.
void MGEdge::compute_continuity2(
	const MGInterval& span,//this edge's binder edge's parameter span.
	int npoint,				//division number of this edge's binder curve interval sspan.
	const MGEdge& edge2,	//the second parameter edge that has the star face.
	double distance[7],		//evaluated data will be set in distance, tangent, and normal.
	double tangent[4],		//See compute_continuity
	double normal[4]
)const{
	MGEdge* be1=binder_edge();
	std::auto_ptr<MGCurve> crv1(be1->curve_limitted());//world curve repr of edge1.

	MGEdge* be2=edge2.binder_edge();
	std::auto_ptr<MGCurve> crv2(be2->curve_limitted());//world curve repr of edge1.

//Here, sx means crv1's parameter value, and tx means crv2's.

	double s_maxd=0., s_mind=0., max_distance=-1., min_distance=-1., mean_distance=0;
	double s_tdif=0., max_tdif=-1., mean_tdif=0.;
	double s_ndif=0., max_ndif=-1., mean_ndif=0.;
	int nspan=npoint-1;
	double rzero2=MGTolerance::rc_zero()*2.;

	MGEdgeMaxDistance max_dist(*crv1,*crv2);
	MGEdgeMinDistance min_dist(*crv1,*crv2);
	MGEdgeTanDifference max_tdiff(*crv1,*crv2); 
	MGEdgeNormalDifference max_ndiff(*this,*crv1,edge2,*crv2); 

	//Find the extreme points in the npoint discrete points.
	double s=span.low_point(), sdelta=span.length().value()/nspan;

	if(sdelta<=crv1->param_span()*rzero2){
		max_dist.set_initial_t(s);
		min_dist.set_initial_t(s);
		max_tdiff.set_initial_t(s);
		max_ndiff.set_initial_t(s);
		max_dist.get_result(distance[0],distance[1],distance[2]);
		min_dist.get_result(distance[3],distance[4],distance[5]);
		distance[6]=distance[2];
		max_tdiff.get_result(tangent[0],tangent[1],tangent[2]);
		tangent[3]=tangent[2];
		max_ndiff.get_result(normal[0],normal[1],normal[2]);
		normal[3]=normal[2];
		return;
	}

	for(int i=0; i<npoint; i++){
		MGPosition P1=crv1->eval(s);
		double t=crv2->closest(P1);

		//max distance info.
		double disti=max_dist.eval(s,t);
		mean_distance+=disti;
		if(max_distance<disti){
			s_maxd=s;
			max_distance=disti;
		}

		//min distance info.
		if(min_distance<0 || min_distance>disti){
			s_mind=s;
			min_distance=disti;
		}

		//max tdif(tangent difference) info.
		double Tanglei=max_tdiff.eval(s,t);
		mean_tdif+=Tanglei;
		if(max_tdif<Tanglei){
			s_tdif=s;
			max_tdif=Tanglei;
		}

		//max ndif(normal difference) info.
		double Nanglei=max_ndiff.eval(s,t);
		mean_ndif+=Nanglei;
		if(max_ndif<Nanglei){
			s_ndif=s;
			max_ndif=Nanglei;
		}

		//increment s to test next point.
		if(i==nspan)
			s=span.high_point();
		else
			s+=sdelta;
	}

	int nrepition;
	double tol=sdelta*rzero2;
	max_dist.solve(s_maxd,sdelta,tol,nrepition);//distance info.
	min_dist.solve(s_mind,sdelta,tol,nrepition);
	max_tdiff.solve(s_tdif,sdelta,tol,nrepition);//tangent difference info.
	max_ndiff.solve(s_ndif,sdelta,tol,nrepition);//normal difference info.

	max_dist.get_result(distance[0],distance[1],distance[2]);
	min_dist.get_result(distance[3],distance[4],distance[5]);
	distance[6]=mean_distance/npoint;
	max_tdiff.get_result(tangent[0],tangent[1],tangent[2]);
	tangent[3]=mean_tdif/npoint;
	max_ndiff.get_result(normal[0],normal[1],normal[2]);
	normal[3]=mean_ndif/npoint;
}
