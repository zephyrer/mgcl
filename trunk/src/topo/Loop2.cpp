/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Interval.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/TrimmedCurve.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/SurfCurve.h"
#include "mg/Surface.h"
#include "mg/CCisect_list.h"
#include "topo/Complex.h"
#include "topo/LPoint.h"
#include "topo/LEPoint.h"
#include "topo/Loop.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/LCisect_vector.h"
#include "topo/LLisect_vector.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGLoop Class.
//MGLoop is a boundary of a face, a boundary of 2D manifold cell.
//MGLoop always accepts parameter space curve and world space curve
//of a boundary curve, and constructs a boundary of a face from the
//two types of curves.

//Input curve direction indicates which part of the face will be target
//part after trimed by the boundary. In 2D space (u,v) of the parameter
//space, LEFT side of the parameter curve along the curve's direction
//is the target part of face.

//Merge loop2 or param_curve to the existing loop data, and build a loop.
//Returned is if merge was done(true) or not(false).
//When no intersection was found, merge is not executed.
//param_curve is parameter space representation of the face
//of which this loop will be a boundary, will make parameter cell.
//world_curve is world coordinate representation of the face
//of which this loop will be a boundary, will make binder cell.
//range1 is parameter range of the curve param_curve.
//When world_curve is input, it will be trimmed according to the range
//of param_curve.
//param_curve and world_curve may be opposite direction.
//When more than one closed loops are detected, first one from
//the old loop start point is employed, and other loops are discarded.
bool MGLoop::merge_trim(const MGCurve& param_curve){
	MGInterval rangep=param_curve.param_range();
	return merge_trim(param_curve,rangep,0);
}
bool MGLoop::merge_trim(const MGCurve& param_curve, const MGInterval& range1){
	return merge_trim(param_curve,range1,0);
}
bool MGLoop::merge_trim(const MGCurve& param_curve, const MGCurve& world_curve){
	assert(star());	//THis must be a member of a cell.

	MGInterval rangep=param_curve.param_range();
	return merge_trim(param_curve,rangep,&world_curve);
}
bool MGLoop::merge_trim(
	const MGCurve& param_curve,
	const MGInterval& range1,
	const MGCurve* world_curve
){
	assert(!world_curve || (world_curve && star()));
	//THis must be a member of a cell when world_curve specified.

	const MGSurface* surf=0;
	if(world_curve)
		surf=surface();

	MGLCisect_vector lcis=
		isect_with_endpoints(MGTrimmedCurve(param_curve, range1));
	//isect of this loop and param_curve.
	//cout<<lcis<<endl;
	size_t n=lcis.entries();
	if(n==0)
		return false;
	if(n>1)
		std::sort(lcis.m_lcisects.begin(),lcis.m_lcisects.end());

	if(closed()){
	//****When trimming closed loop.
	bool in=is_inner_boundary();

	if(n<2)
		return false;
	double ang1,ang2;
	MGVector vec1, vec2;
	MGLCisect lci1=lcis.m_lcisects[0], lci2=lcis.m_lcisects[1];
	vec1=eval(lci1.lp(),1);
	vec2=param_curve.eval(lci1.t(),1);
	ang1=(vec1*vec2).ref(2);

	vec1=eval(lci2.lp(),1);
	vec2=param_curve.eval(lci2.t(),1);
	ang2=(vec1*vec2).ref(2);
	if(ang1*ang2>0.) return false;

	MGLCisect lcisave;
	double s1,s2;
	if(ang1<0.){
		s1=lci1.t(); s2=lci2.t();
	}else{
		s1=lci2.t(); s2=lci1.t();
		lcisave=lci1; lci1=lci2; lci2=lcisave;
	}
	//(lci1,s1)=parameter of ang<0. and (lci2,s2)= of ang>0.

	if(!in && s1<s2)
		return false;///Case tha trimming out of outer loop.

	trim(lci1,lci2);
	MGEdge* add_edge;
	if(in && s1<s2){
	//Case that trimming out of inner loop and becomming open loop.
		MGInterval rng1(param_curve.param_s(),s1);
		if(surf)
			add_edge=new MGEdge(*surf,param_curve,rng1,*world_curve);
		else
			add_edge=new MGEdge(param_curve,rng1);
		first_edge()->join(true,add_edge);
		MGInterval rng2(s2,param_curve.param_e());
		if(surf)
			add_edge=new MGEdge(*surf,param_curve,rng2,*world_curve);
		else
			add_edge=new MGEdge(param_curve,rng2);
		last_edge()->join(false,add_edge);//	cout<<(*this)<<endl;
	}else{
	//Case that becomming closed loop.
		MGInterval rng(s2,s1);
		if(surf) add_edge=new MGEdge(*surf,param_curve,rng,*world_curve);
		else	 add_edge=new MGEdge(param_curve,rng);
		last_edge()->join(false,add_edge);//	cout<<(*this)<<endl;
		last_edge()->join(false,first_edge());//	cout<<(*this)<<endl;
	}
	return true;

	}else{
	//****When trimming open loop.

	MGLoop save;
	size_t i=0;
	MGLCisect lci=lcis.m_lcisects[0], lci2;
	double t1, t2=lci.t();
	MGLEPoint let1(lci);
	if(!let1.is_end_point()){

	MGVector vec1=eval(lci.lp(),1), vec2=param_curve.eval(t2,1);
	if(let1.is_start_point() || (vec1*vec2).ref(2)<0.){
		//start of loop and end of param_curve.
		if(n>1){
			lci2=lcis.m_lcisects[1];
			t1=lci2.t();
		}
		if(n>1 && t1<t2){
		//Loop makes a closed loop.
			trim(lci,lci2);
			MGEdge* add_edge;
			MGInterval rng(t1,t2);
			if(surf)
				add_edge=new MGEdge(*surf,param_curve,rng,*world_curve);
			else
				add_edge=new MGEdge(param_curve,rng);
			first_edge()->join(true,add_edge);
			last_edge()->join(false,add_edge);
			return true;
		}else{
			if(n>1){
				save=*this;
				trim(lci,lci2);
			}else trim_start(lci);
			MGInterval rng(range1);
			rng.set_high_point(t2);
			MGEdge* add_edge;
			if(surf)
				add_edge=new MGEdge(*surf,param_curve,rng,*world_curve);
			else
				add_edge=new MGEdge(param_curve,rng);
			first_edge()->join(true,add_edge);
			if(n==1) return true;
			i=1;
		}
	}else{
		if(n>1) save=*this;
		trim_end(lci);
	}

	}

	MGLoop span, *spanp;
	MGEReal rend=range1.high();
	while(n>i){
		i++;
		t2=lci.t();
		MGEReal e1(t2),e2(rend);
		if(n>i){
			lci=lcis.m_lcisects[i];
			e2=MGEReal(lci.t());
		}
		MGInterval rng(e1,e2);
		MGEdge* add_edge;
		if(surf)
			add_edge=new MGEdge(*surf,param_curve,rng,*world_curve);
		else
			add_edge=new MGEdge(param_curve,rng);
		last_edge()->join(false,add_edge);
		if(n>i){
			i++;
			if(n>i){
				lci2=lcis.m_lcisects[i];
				span=save; span.trim(lci,lci2);
				spanp=&span;
				lci=lci2;
			}else{
				save.trim_start(lci);
				spanp=&save;
			}
			join(false,*spanp);
		}
	}
	}
	return true;
}

bool MGLoop::merge_trim(const MGLoop& loop2){
	bool merged;
	if(closed())
		merged=merge_loop(loop2);
	else if(loop2.closed()){
		MGLoop temp(loop2);
		if(merged=temp.merge_loop(*this))
			*this=temp;
	}else
		merged=merge_loop(loop2);
	return merged;
}
bool MGLoop::merge_loop(const MGLoop& loop2){
	MGLLisect_vector llis=isect(loop2);	//isect of this loop and loop2.

	size_t n=llis.entries();
	if(n==0) return false;

	std::sort(llis.m_llisects.begin(),llis.m_llisects.end());
if(closed()){

//****When trimming closed loop.
	if(n<2)
		return false;

	bool in=is_inner_boundary();

	double ang1,ang2;
	MGVector vec1, vec2;
	MGLLisect lli1=llis.m_llisects[0], lli2=llis.m_llisects[1];
	vec1=eval(lli1.isect1(),1);
	vec2=loop2.eval(lli1.isect2(),1);
	ang1=(vec1*vec2).ref(2);
	vec1=eval(lli2.isect1(),1);
	vec2=loop2.eval(lli2.isect2(),1);
	ang2=(vec1*vec2).ref(2);
	if(ang1*ang2>0.) return false;

	MGLLisect llisave;
	if(ang2<0.){ llisave=lli1; lli1=lli2; lli2=llisave;}
	//lli1=parameter of ang<0. and lli2= of ang>0.

	MGLoop loop2_t=loop2;
	if(loop2.closed()){
		loop2_t.trim(lli2.isect2(),lli1.isect2());
		trim(lli1.isect1(), lli2.isect1());
		join(true,loop2_t);
		first_edge()->join(true,last_edge());
	}else{
		if(!in && lli1.isect2()<lli2.isect2()) return false;
			///Case tha trimming out of outer loop.

		trim(lli1.isect1(),lli2.isect1());
		if(in && lli1.isect2()<lli2.isect2()){
		//Case that trimming out of inner loop and becomming open loop.
			MGLoop loop2_t2=loop2;
			loop2_t2.trim_end(lli1.isect2());
			join(true,loop2_t2);
			loop2_t.trim_start(lli2.isect2());
			join(false,loop2_t);
		}else{
		//Loop makes a closed loop.
			loop2_t.trim(lli2.isect2(),lli1.isect2());
			join(false,loop2_t);				//cout<<(*this)<<endl;
			first_edge()->join(true,last_edge());
		}
	}
	return true;

}else{
//****When trimming open loop.

	MGLoop save;	//To save *this(if necessary).
	size_t i=0;
	MGLLisect lli=llis.m_llisects[0], lli2;
	MGLEPoint t1=lli.isect1(), t2=lli.isect2();
	if(!t1.is_end_point()){

	MGVector vec1=eval(t1,1), vec2=loop2.eval(t2,1);
	if(t1.is_start_point() || (vec1*vec2).ref(2)<0.){
		//Yes, when lower part of loop2 is on the left side of this loop.
		if(n>1){
			lli2=llis.m_llisects[1];
			t1=lli2.isect2();
		}
		if(n>1 && t1<t2){
		//Loop makes a closed loop.
			trim(lli.isect1(),lli2.isect1());
			MGLoop loop2_t=loop2;loop2_t.trim(t1,t2);
			join(false,loop2_t);
			first_edge()->join(true,last_edge());
			return true;
		}else{
			if(n>1){
				save=*this;
				trim(lli.isect1(),lli2.isect1());
			}else trim_start(lli.isect1());

			MGLoop loop2_t=loop2;loop2_t.trim_end(lli.isect2());
			join(true,loop2_t);
			if(n==1) return true;
			lli=lli2; i=1;
		}
	}else{
		//Yes, when higher part of loop2 is on the left side of this loop.
		if(n>1) save=*this;
		trim_end(lli.isect1());
	}

	}

	MGLoop span, *spanp;
	while(n>i){
		i++;
		MGLoop loop2_t=loop2;
		if(n>i){
			lli2=llis.m_llisects[i];
			loop2_t.trim(lli.isect2(),lli2.isect2());
			lli=lli2;
		}else{
			loop2_t.trim_start(lli.isect2());
		}
		join(false,loop2_t);
		if(n>i){
			i++;
			if(n>i){
				lli2=llis.m_llisects[i];
				span=save; span.trim(lli.isect1(),lli2.isect1());
				spanp=&span;
				lli=lli2;
			}else{
				save.trim_start(lli.isect1());
				spanp=&save;
			}
			join(false,*spanp);
		}
	}
}
	return true;
}

//Trim the loop. Result loop is from t1 to t2;
//The loop can be closed one. In this case, t1 can be >t2.
//When not closed, t1 must be less than t2.
void MGLoop::trim(const MGLEPoint& t1, const MGLEPoint& t2){
	assert(closed()||(!closed()&&t1<t2));
		//When not closed loop, t1 must be <t2.

	if(t1<t2){
		trim_end(t2);
		trim_start(t1);
	}else{
	//When closed loop and trimming across connecting point.
		MGLoop lp2=*this;
		MGLEPoint t12(*this,t1,lp2);//cout<<t12;
		trim_end(t2);//cout<<*this;//**********
		lp2.trim_start(t12);//cout<<lp2;//***********
		join(true,lp2);//cout<<*this<<endl;//******
	}
}

//Trim the loop. Result is from start to t2.
void MGLoop::trim_end(const MGLEPoint& t2){
	size_t n=number_of_pcells();
	for(size_t i=t2.edge_num()+1; i<n;i++)
		erase_last_pcell();
	trim_param_set(t2,false);
}

//Dedicated function of trim_start, end, will set new parameter
//range of the trimming edge;
//When start=true, set as start point data, else end point data.
void MGLoop::trim_param_set(const MGLEPoint& t, bool start){
	MGEdge* e=t.edge_to_update();
	e->trim(t.param(),start);
	m_area=0.;
	m_box.set_null();
}

//Trim the loop. Result is from t1 to end.
void MGLoop::trim_start(const MGLEPoint& t1){
	size_t n=number_of_pcells();
	size_t en=t1.edge_num();
	for(size_t i=0; i<en; i++)
		erase_first_pcell();
	trim_param_set(t1,true);
}
