/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <iomanip>
#include "mg/Interval.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/Tolerance.h"
#include "topo/LEPoint.h"
#include "topo/Loop.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "topo/LCisect_vector.h"
#include "topo/LLisect_vector.h"

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

struct LCComp{
	bool operator()(const MGLCisect& lc1, const MGLCisect& lc2){return lc1.t()<lc2.t();};
};

//Merge with param_curve as a network loop.
//Function's return value is:
//true: merge is done and param_curve is processed into this loop as a network.
//false: merge is not performed since no intersection with this loop were found.
//When false is returned, this loop does not have the input param_curve information
//as an edge information.
//This loop must not be empty loop, and the kind is always changed to NETWORK.
bool MGLoop::merge_network(const MGCurve& param_curve){
	m_kind=NETWORK;

	MGLCisect_vector lcis=isect_with_endpoints(param_curve);
	//isect of this loop and param_curve.
	size_t n=lcis.entries();
	if(n==0)
		return false;

	if(n>1)
		std::sort(lcis.m_lcisects.begin(),lcis.m_lcisects.end(), LCComp());
	//std::cout<<lcis<<std::endl;/////////////****************

	MGEdge* pre=0; MGEdge* aft=0;
	const MGLCisect& lc_start=lcis[0];
	double tnow=param_curve.param_s(), tspan=param_curve.param_span(),
		tnext=lc_start.t();

	if(!MGREqual_base(tnext,tnow,tspan)){
		const MGLEPoint& le=lc_start.lp();
		if(make_vertex(le,pre,aft))
			lcis.update_lepoint(le);
		//std::cout<<(*this)<<std::endl;//************************
		MGEdge* e=new MGEdge(param_curve, MGInterval(tnow,tnext));
		if(pre)
			pre->connect_at_id(1,e,1);
		else
			aft->connect_at_id(0,e,1);
		//std::cout<<(*this)<<std::endl;//************************
		//if(n<=1)//Case of only one intersection point.
		//	return true;
		tnow=tnext;
	}

	size_t nm1=n-1;
	for(size_t i=0; i<nm1; i++){
		//Process one edge of param_curve(the start and end point connect).
		const MGLCisect& lc_pre=lcis[i];
		if(!pre){
			const MGLEPoint& le=lc_pre.lp();
			if(make_vertex(le,pre,aft))
				lcis.update_lepoint(le);
		}
		const MGLCisect& lc_aft=lcis[i+1];

		tnext=lc_aft.t();
		MGEdge* ei1=new MGEdge(param_curve, MGInterval(tnow,tnext));		

		if(pre)
			pre->connect_at_id(1,ei1,0);
		else
			aft->connect_at_id(0,ei1,0);

		const MGLEPoint& le=lc_aft.lp();
		if(make_vertex(le,pre,aft))
			lcis.update_lepoint(le);
		if(pre)
			pre->connect_at_id(1,ei1,1);
		else
			aft->connect_at_id(0,ei1,1);
		tnow=tnext;
		//std::cout<<(*this)<<std::endl;//************************
	}

	const MGLCisect& lc_end=lcis[nm1];
	tnext=param_curve.param_e();
	if(!MGREqual_base(tnext,tnow,tspan)){
		if(!pre){
			const MGLEPoint& le=lc_end.lp();
			if(make_vertex(le,pre,aft))
				lcis.update_lepoint(le);
		}
		MGEdge* e=new MGEdge(param_curve, MGInterval(tnow,tnext));
		if(pre)
			pre->connect_at_id(1,e,0);
		else
			aft->connect_at_id(0,e,0);
	}
	//std::cout<<(*this)<<std::endl;//************************
	return true;
}

//Remove pendent edge.
void MGLoop::remove_pendent_edge(const MGFSurface& face){
	pcellItr i=pcell_begin(), iend=pcell_end();
	for(;i!=iend;){
		MGEdge* edge=static_cast<MGEdge*>(*i);i++;
		MGEdge* pre=edge->pre_edge();
		MGEdge* aft=edge->aft_edge();
		if(pre&&aft)
			continue;
		if(!pre){
			MGPosition uv=edge->start_point();
			int in=face.in_range_with_on(uv);
			if(0<=in && in<=2){
				erase_pcell(edge);
				continue;
			}
		}
		if(!aft){
			MGPosition uv=edge->end_point();
			int in=face.in_range_with_on(uv);
			if(0<=in && in<=2){
				erase_pcell(edge);
				continue;
			}
		}
	}
}

class mgCPoints;
class mgCPointsVec;
class mgCONNECT_POINT{
public:
	enum OUTIN{coming_in, going_out};
private:
	MGTrimLoop* m_trim_loop;//MGTrimLoop that touch the inner loop.
	OUTIN m_outin;
		//indicates if the contact point is coming_in of m_trim_loop(that is,
		//the end point of m_trim_loop), or going_out(that is, the start point
		//of m_trim_loop).

public:
	mgCONNECT_POINT():m_trim_loop(0){;};
	mgCONNECT_POINT(MGTrimLoop* loop, OUTIN outin)
		:m_trim_loop(loop), m_outin(outin){;}

	bool operator< (const mgCONNECT_POINT& linf2)const;
	bool operator== (const mgCONNECT_POINT& linf2)const;
	bool operator!= (const mgCONNECT_POINT& linf2)const{return !operator==(linf2);};
	bool is_null()const{
		return m_trim_loop->is_null();
	}
	bool is_going_out(){
		return m_outin==going_out;
	};
	bool is_coming_in(){
		return m_outin==coming_in;
	};
	MGLEPoint lep(){
		if(m_outin==coming_in)
			return m_trim_loop->end_lep();
		else
			return m_trim_loop->start_lep();
	};
	std::auto_ptr<MGLoop> loop_clone(){
		return std::auto_ptr<MGLoop>(new MGLoop(*(m_trim_loop->loop())));
	};
	OUTIN outin()const{return m_outin;};
	MGTrimLoop* trim_loop(){return m_trim_loop;};
	const MGTrimLoop* trim_loop()const{return m_trim_loop;};

	//obtain the next point of this connect point.
	//The next point may be in the same loop, or in a different loop.
	//function's return value is the mgCONNECT_POINT of the next point.
	mgCONNECT_POINT next_point(
		mgCPointsVec& cpointsVec//container of MGTrimLoop.
	);

	//obtain the previous point of this connect point.
	//The next point may be in the same loop, or in a different loop.
	//function's return value is the mgCONNECT_POINT of the previous point.
	mgCONNECT_POINT prev_point(
		mgCPointsVec& cpointsVec//container of MGTrimLoop.
	);

	//Debug Function
	friend std::ostream& operator<< (std::ostream& out, const mgCONNECT_POINT& cpoint);

};

//Debug Function
std::ostream& operator<< (std::ostream& out, const mgCONNECT_POINT& cpoint){
	out<<"mgCONNECT_POINT::m_trim_loop="<<*(cpoint.m_trim_loop);
	if(cpoint.m_outin==mgCONNECT_POINT::coming_in){
		out<<",coming_in";
	}else{
		out<<",going_out";
	}
	return out;
}

class mgCPoints{
public:
	std::vector<mgCONNECT_POINT> m_cpoints;
	size_t size()const{return m_cpoints.size();}
	void sort(){std::sort(m_cpoints.begin(), m_cpoints.end());}
	void push_back(const mgCONNECT_POINT& cpoint){m_cpoints.push_back(cpoint);};

	int next_id(int j)const{return (j+1)%m_cpoints.size();};
	int prev_id(int j)const{
		int n=m_cpoints.size();
		return (j+n-1)%n;
	};
	mgCONNECT_POINT next_point(int j){
		return m_cpoints[next_id(j)];
	}
	mgCONNECT_POINT prev_point(int j){
		return m_cpoints[prev_id(j)];
	}

	//get the id of m_cpoints that includes tloop
	int find_trim_loop(const MGTrimLoop* tloop, mgCONNECT_POINT::OUTIN out_in)const;

	//Debug Function
	friend std::ostream& operator<< (std::ostream& out, const mgCPoints& cpoints);
};

class mgCPointsVec{

const MGFace& m_face;//Target face to trim
const MGPosition& m_uv;//trimming positioning data.
MGPvector<MGTrimLoop> m_trim_loops;//Trimloop's
std::vector<mgCPoints> m_cpoints_vec;//vector of mgCPoints.

public:

explicit mgCPointsVec(
	const MGFace& face, //target face.
	const MGPosition& uv=mgNULL_VEC
):m_face(face), m_uv(uv), m_cpoints_vec(face.get_number_of_boundaries()){;};

//Extract MGTrimLoop & mgCPoints info, and build this mgCPointsVec data.
//While processing, a closed outer loop or closed inner loop is found, store it in
//the oloops or iloops.
//When m_uv is not null, oloops contains the only closed outer loop that contain m_uv,
//and return function's return code as 1. If the closed loop was not found,
//generate m_trim_loops(the vector of MGTrimLoop).
//When m_uv is null, put all the outer loops into oloops, and inner loops into iloops.
//Function's return value is
//  1: when m_uv is not null, a closed loop that includes m_uv was extracted into oloops.
//  0: when m_uv is not null, no closed outer loop that includes m_uv was found
//     and trim loop info were output into m_trim_loops and
//     mgCONNECT_POINT's into m_cpoints_vec. Closed inner loops may be output to iloops.
//When m_uv is null, function's return value is always 0.
int extract_trim_loop_info(
	const MGLoop& network,	//original loops to trim.
	MGPvector<MGLoop>& oloops,//When m_uv is not nul, the closed outerboundary loop
		//that includes m_uv is output and the return code ==1.
		//When m_uv is null, all the detected closed outerboundary loops will be output,
		//and the function's return value is 0.
	MGPvector<MGLoop>& iloops//closed inner boundary loop that
							//includes uv is output when return code =0 or 1.
);

//Extract a boundary out of networks that includes m_uv.
//m_uv is a parameter of m_face.
//Function's return value is
//  0: no loops were extracted,
//  2: Loops to trim face were output into used_tloops.
//     If trimming is performed by used_tloops[i] for i=0,...,n-1 one by one,
//     face will be trimmed by the part of the face that includes uv.
int extract_uv_boundary(
	std::auto_ptr<MGLoop>& loop,	//the loop that includes m_uv will be output
		//when return code=2.
	std::deque<MGTrimLoop*>& used_tloops
);

//extract the loop whose start point is first_point.
//Function's return value is:
//0: no loop was extracted since first_point was null.
//1: a loop was extracted into loop.
int extract_loop(
	mgCONNECT_POINT& first_point,
	std::auto_ptr<MGLoop>& loop, //the loop will be output.
	std::deque<MGTrimLoop*>& tloops//used MGTrimLoop to extract loop will be output.
);

mgCPoints& operator[](size_t i){return m_cpoints_vec[i];};
void push_at(size_t i, const mgCONNECT_POINT& cpoint){m_cpoints_vec[i].push_back(cpoint);};
size_t size()const{return m_cpoints_vec.size();};
bool no_uv()const{return m_uv.is_null();};

friend std::ostream& operator<< (std::ostream& out, const mgCPointsVec& cpvec);

};

std::ostream& operator<< (std::ostream& out, const mgCPointsVec& cpvec){
	out<<"mgCPointsVec::"<<std::endl;
	size_t n=cpvec.m_trim_loops.size();
	for(size_t i=0; i<n; i++)
		out<<" m_trim_loops["<<i<<"]:"<<*(cpvec.m_trim_loops[i]);
	return out;
}

bool mgCONNECT_POINT::operator< (const mgCONNECT_POINT& cpoint2)const{
	MGLEPoint le1;
	if(m_outin==coming_in)
		le1=m_trim_loop->end_lep();
	else
		le1=m_trim_loop->start_lep();

	MGLEPoint le2;
	if(cpoint2.m_outin==coming_in)
		le2=cpoint2.m_trim_loop->end_lep();
	else
		le2=cpoint2.m_trim_loop->start_lep();

	if(le1 == le2){
		if(m_outin==going_out)
			return true;
		else
			return false;
	}
	return le1<le2;
}
bool mgCONNECT_POINT::operator== (const mgCONNECT_POINT& cpoint2)const{
	if(m_trim_loop!=cpoint2.m_trim_loop)
		return false;
	if(m_outin!=cpoint2.m_outin)
		return false;
	return true;
}

//obtain the next point of this connect point.
//The next point may be in the same loop, or in a different loop.
//function's return value is the mgCONNECT_POINT of the next point.
mgCONNECT_POINT mgCONNECT_POINT::next_point(
	mgCPointsVec& cpointsVec
){
	MGTrimLoop* tloop=trim_loop();
	size_t loopid=tloop->end_loopid();
	mgCPoints& cpoints=cpointsVec[loopid];
	int id=cpoints.find_trim_loop(tloop,coming_in);
	if(is_coming_in()){
		//The next point of the end point of tloop.
		return cpoints.next_point(id);
	}else{
		//The next point of the start point of trim_loop
		//(that is the end point of tloop).
		return cpoints.m_cpoints[id];
	}
}

//obtain the previous point of this connect point.
//The previous point may be in the same loop, or in a different loop.
//function's return value is the id of cpoints of the previous point.
mgCONNECT_POINT mgCONNECT_POINT::prev_point(
	mgCPointsVec& cpointsVec
){
	MGTrimLoop* tloop=trim_loop();
	size_t loopid=tloop->start_loopid();
	mgCPoints& cpoints=cpointsVec[loopid];
	int id=cpoints.find_trim_loop(tloop,going_out);
	if(is_coming_in()){
		//The previous point of the end point of trim_loop
		//(that is the start point of tloop).
		return cpoints.m_cpoints[id];
	}else{
		//The previous point of the start point of tloop.
		return cpoints.prev_point(id);
	}
}

//extract the loop whose start point is first_point.
//Function's return value is:
//0: no loop was extracted since first_point was null.
//1: a loop was extracted into loop.
int mgCPointsVec::extract_loop(
	mgCONNECT_POINT& first_point,
	std::auto_ptr<MGLoop>& loop, //the loop will be output.
	std::deque<MGTrimLoop*>& tloops//used MGTrimLoop to extract loop will be output.
){
	if(first_point.is_null())
		return 0;

	mgCONNECT_POINT cpoint, npoint;
	if(first_point.is_coming_in()){
		cpoint=first_point;
		loop=std::auto_ptr<MGLoop>(new MGLoop);
	}else{		
		loop=first_point.loop_clone();
		tloops.push_back(first_point.trim_loop());
		cpoint=first_point.next_point(*this);
	}

	//construct loop to the forward direction.
	size_t ntloops=m_trim_loops.size();
	for(size_t counter=0; counter<ntloops; counter++){//This counter is to avoid an infinite loop.
		assert(cpoint.is_coming_in());
		npoint=cpoint.next_point(*this);assert(npoint.is_going_out());

		if(npoint.is_null())
			return 0;//This must not happen.

		MGLEPoint ts=cpoint.lep(), te=npoint.lep();
		assert(ts.loop()==te.loop());
		std::auto_ptr<MGLoop> loop2=trim_out_subloop(ts,te);
		loop->join(false,loop2);
		if(npoint==first_point){
			loop->make_close();
			return 1;
		}

		loop->join(false,npoint.loop_clone());
		tloops.push_back(npoint.trim_loop());
		cpoint=npoint.next_point(*this);
		if(cpoint==first_point){
			loop->make_close();
			return 1;
		}
	}
	return 0;//This must not happen.
}

//Build a boundary loop by tracking tloops backward.
//Function's return value is:
//0: no loops were extracted
//1: a boundary loop was extracted from the first non-null tloops
int build_boundary_loop(
	std::deque<MGTrimLoop*>& tloops,//used MGTrimLoop to extract loop is input.
	std::auto_ptr<MGLoop>& loop,//the loop will be output.
	std::vector<bool>& used_loops//vector of length face.number_of_boundaries().
		//used loops to build loop will be set true.
){
	if(tloops.empty())
		return 0;
	MGTrimLoop* tloop;
	while(!tloops.empty()){
		tloop=tloops.back(); tloops.pop_back();
		if(!tloop->is_null())
			break;
	}
	if(tloop->is_null())
		return 0;

	size_t nboundaries=used_loops.size();
	for(size_t j=0; j<nboundaries; j++)
		used_loops[j]=false;
	tloop->set_used_loop_flag(used_loops);

	size_t first_point_loopid=tloop->end_loopid();
	size_t next_point_loopid=tloop->start_loopid();
	MGLEPoint first_point=tloop->end_lep();
	MGLEPoint ts=tloop->start_lep();
	loop=std::auto_ptr<MGLoop>(tloop->release_loop());
	loop->negate();
	while(next_point_loopid!=first_point_loopid && !tloops.empty()){
		//while next_point_loopid!=first_point_loopid, tloops must not be empty.
		tloop=tloops.back(); tloops.pop_back();
		MGLEPoint te=tloop->end_lep();
		assert(ts.loop()==te.loop());
		std::auto_ptr<MGLoop> loop2=trim_out_subloop(ts,te);
		loop->join(false,loop2);

		next_point_loopid=tloop->start_loopid();used_loops[next_point_loopid]=true;
		ts=tloop->start_lep();
		MGLoop* loop3=tloop->release_loop();
		loop3->negate();
		loop->join(false,loop3);
	}
	std::auto_ptr<MGLoop> loop4=trim_out_subloop(ts,first_point);
	loop->join(false,loop4);
	loop->make_close();

	return 1;
}

//get the id of m_cpoints that includes tloop
int mgCPoints::find_trim_loop(const MGTrimLoop* tloop, mgCONNECT_POINT::OUTIN out_in)const{
	int n=m_cpoints.size();
	for(int i=0; i<n; i++){
		const mgCONNECT_POINT& cpi=m_cpoints[i];
		if(tloop==cpi.trim_loop())
			if(cpi.outin()==out_in)
				return i;
	}
	return -1;
}

//Debug Function
std::ostream& operator<< (std::ostream& out, const mgCPoints& cpoints){
	int n=cpoints.m_cpoints.size();
	out<<"mgCPoints, size="<<n<<",";
	for(int i=0; i<n; i++){
		out<<std::endl<<i<<"::"<<cpoints.m_cpoints[i]<<std::endl;
	}
	return out;
}

//Extract MGTrimLoop & mgCPoints info, and build this mgCPointsVec data.
//While processing, a closed outer loop or closed inner loop found will be stored in
//the oloops or iloops.
//When m_uv is not null, oloops contains the only closed outer loop that contain m_uv,
//and return function's return code as 1. If the closed loop was not found,
//generate m_trim_loops(the vector of MGTrimLoop).
//When m_uv is null, put all the outer loops into oloops, and inner loops into iloops.
//Function's return value is
//  1: when m_uv is not null, a closed loop that includes m_uv was extracted into oloops.
//  0: when m_uv is not null, no closed outer loop that includes m_uv was found
//     and trim loop info were output into m_trim_loops and
//     mgCONNECT_POINT's into m_cpoints_vec. Closed inner loops may be output to iloops.
//When m_uv is null, function's return value is always 0.
int mgCPointsVec::extract_trim_loop_info(
	const MGLoop& network,	//original loops to trim.
	MGPvector<MGLoop>& oloops,//When m_uv is not nul, the closed outerboundary loop
		//that includes m_uv is output and the return code ==1.
		//When m_uv is null, all the detected closed outerboundary loops will be output,
		//and the function's return value is 0.
	MGPvector<MGLoop>& iloops//closed inner boundary loop that
							//includes uv is output when return code =0 or 1.
){
	size_t nedge=network.number_of_edges();
	if(!nedge)
		return 0;//Any loop was not found.

	enum SFLAG{
		NOT_YET=0,
		SAME=1,
		OPPOSITE=2,
		BOTH=3
	};
	std::vector<int> searched(nedge,NOT_YET);
		//NOT_YET: not searched yet for the network.edge(i),
		//SAME: searched for the same direction of the original edge,
		//OPPOSITE: searched for the opposite direction of the original
		//BOTH(=SAME+OPPOSITE): searched for the both directions.

	//Serach for the same direction.
	bool handled=true;
	while(handled){
	handled=false;
	for(size_t i=0; i<nedge; i++){
		int ei_flag=searched[i];
		if(ei_flag>=BOTH)//If searched in the both directions.
			continue;

		std::auto_ptr<MGLoop> loop2(new MGLoop);
		handled=true;
		const MGEdge* first_edge=network.edge(i);//save the 1st edge to judge loop to the 1st.
		bool first_was_same_direction;
			//True means the target loop has the same direction as the previous edge.

		if(ei_flag==OPPOSITE || ei_flag==NOT_YET){//If search was not done at the same
			first_was_same_direction=true;
		}else{//If search was not done at the opposite
			first_was_same_direction=false;
		}

		const MGEdge* edgei=first_edge;
		bool is_same_direction=first_was_same_direction;
		size_t edge_id=i;
		const MGEdge* next_edge=0;
		do{//Search is performed toward the loop direction.
			if(is_same_direction){
				searched[edge_id]+=SAME;
				loop2->append(edgei->clone());
				next_edge=edgei->aft_edge();
			}else{
				searched[edge_id]+=OPPOSITE;
				MGEdge* ei=edgei->clone();
				ei->negate();
				loop2->append(ei);
				next_edge=edgei->aft_edge(false);//aft edge at start point of the edge.
			}

			if(next_edge && next_edge!=first_edge){
				if(is_same_direction){
					if(edgei->is_connected_and_same_direction(false,*next_edge))
						is_same_direction=true;
					else
						is_same_direction=false;
				}else{
					if(edgei->is_connected_and_same_direction(true,*next_edge))
						is_same_direction=false;
					else
						is_same_direction=true;
				}
				edgei=next_edge;
				edge_id=network.edge_num(edgei);
			}
		}while(next_edge!=0 && next_edge!=first_edge);

		if(next_edge==first_edge){
			//std::cout<<*loop2<<std::endl;//////********
			loop2->make_close();
			if(no_uv()){
				if(loop2->is_outer_boundary()){
					oloops.push_back(loop2.release());
				}else
					iloops.push_back(loop2.release());
			}else if(loop2->inside(m_uv)){
				if(loop2->is_outer_boundary()){
					oloops.push_back(loop2.release());
					return 1;
				}else
					iloops.push_back(loop2.release());
			}
		}else{
			assert(!next_edge);
			//When edge reached to the boundary of the original face,
			//opposite direction searching must be done.
			if(first_was_same_direction){
				next_edge=first_edge->pre_edge();
			}else{
				next_edge=first_edge->pre_edge(false);//pre edge at the end point of 1st edge.
			}
			edgei=first_edge;
			is_same_direction=first_was_same_direction;
			while(next_edge){//Search is performed against the loop direction.
				if(is_same_direction){
					if(edgei->is_connected_and_same_direction(true,*next_edge))
						is_same_direction=true;
					else
						is_same_direction=false;
				}else{
					if(edgei->is_connected_and_same_direction(false,*next_edge))
						is_same_direction=false;
					else
						is_same_direction=true;
				}
				edgei=next_edge;
				edge_id=network.edge_num(edgei);
				if(is_same_direction){
					searched[edge_id]+=SAME;
					loop2->prepend(edgei->clone());
					next_edge=edgei->pre_edge();
				}else{
					searched[edge_id]+=OPPOSITE;
					MGEdge* ei=edgei->clone();
					ei->negate();
					loop2->prepend(ei);
					next_edge=edgei->pre_edge(false);//pre edge at the end point of the edge.
				}
			}
			MGPosition uvS=loop2->start_point();
			MGPosition uvE=loop2->end_point();
			int inS=m_face.in_range_with_on(uvS),
				inE=m_face.in_range_with_on(uvE);
			int lidS=0; if(inS<0) lidS=-inS;
			int lidE=0; if(inE<0) lidE=-inE;
			double d;
			MGLEPoint leS=m_face.loop(size_t(lidS))->closest(uvS,d),
					  leE=m_face.loop(size_t(lidE))->closest(uvE,d);
			//std::cout<<leS<<","<<leE<<std::endl;//////////*********
			MGLoop* lp=loop2.release();
			MGTrimLoop* trloop=new MGTrimLoop(lp,inS,leS,inE,leE);
			m_trim_loops.push_back(trloop);
			push_at(lidS,mgCONNECT_POINT(trloop,mgCONNECT_POINT::going_out));
			push_at(lidE,mgCONNECT_POINT(trloop,mgCONNECT_POINT::coming_in));
		}

	}
	}

	size_t nboundaries=size();
	for(size_t i=0; i<nboundaries; i++){
		mgCPoints& lveci=(*this)[i];
		lveci.sort();
	}

	return 0;
}

//Extract a boundary out of networks that includes m_uv.
//m_uv is a parameter of m_face.
//Function's return value is
//  0: no loops were extracted,
//  2: Loops to trim face were output into used_tloops.
//     If trimming is performed by used_tloops[i] for i=0,...,n-1 one by one,
//     face will be trimmed by the part of the face that includes uv.
int mgCPointsVec::extract_uv_boundary(
	std::auto_ptr<MGLoop>& loop,	//the loop that includes m_uv will be output
		//when return code=2.
	std::deque<MGTrimLoop*>& used_tloops
){
	size_t nboundaries=size();
	for(size_t i=0; i<nboundaries; i++){
		mgCPoints& lveci=(*this)[i];
		size_t nlv=lveci.size();
		assert((nlv%2)==0);//nlv must be even(pair of coming_in and going_out).
		if(!nlv)
			continue;

		for(size_t j=0; j<nlv; j++){
			mgCONNECT_POINT& fpoint=lveci.m_cpoints[j];
			int extract_inf=extract_loop(fpoint,loop,used_tloops);
			if(!extract_inf)
				continue;

			//std::cout<<(*loop)<<std::endl;////************
			size_t ntloop=used_tloops.size(), k;
			if(!loop->inside(m_uv)){
				for(k=0; k<ntloop; k++)
					used_tloops[k]->set_null();
				continue;//if loop does not includes uv.
			}
			return 2;
		}
	}
	return 0;
}

//Trim the face giving networks loops. Trimming is done by removing the smallest
//closed area out of networks that includes the parameter position uv(u,v).
void MGFace::trim(
	const MGPvector<MGLoop>& networks,
	const MGPosition& uv,
	MGPvector<MGFace>& faces//Result trimmed face(s) will be appended.
)const{
	mgCPointsVec cpointsVec(*this,uv);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	MGPvector<MGLoop> iloops;
	int nnetworks=networks.size();
	for(int i=nnetworks-1; i>=0; i--){
		MGPvector<MGLoop> oloops;
		int closed_loop_obtained
			=cpointsVec.extract_trim_loop_info(*(networks[i]),oloops,iloops);
		if(closed_loop_obtained){
			//std::cout<<" Closed loop was extracted."<<std::endl;
			//std::cout<<(*cloop)<<std::endl;
			oloops[0]->negate();
			MGFace* face2=new MGFace(*this);
			face2->add_boundary(oloops.release(0));
			faces.push_back(face2);
			return;
		}
	}

	std::auto_ptr<MGLoop> loop;
	std::deque<MGTrimLoop*> used_tloops;
	int eout=cpointsVec.extract_uv_boundary(loop,used_tloops);
	if(!eout){
		size_t niloops=iloops.size();
		for(size_t i=0; i<niloops; i++){
			MGLoop* iloop=iloops.release(i);
			iloop->negate();
			MGFace* face2=new MGFace(*this);
			face2->add_boundary(iloop);
			faces.push_back(face2);
		}
		return;
	}

	/*std::cout<<" Trim_loops were extracted."<<std::endl;
	int ntloops=used_tloops.size();
	std::cout<<"Num of trim loops="<<ntloops<<std::endl;
	for(int i=0; i<ntloops; i++){
		MGTrimLoop& tl=*(used_tloops[i]);
		if(!tl.is_null())
			std::cout<<","<<*(used_tloops[i]);
	}
	std::cout<<std::endl;*/

	size_t nboundaries=number_of_boundaries();
	std::vector<bool> used_loops(nboundaries);
	std::auto_ptr<MGLoop> new_loop;
	while(build_boundary_loop(used_tloops,new_loop,used_loops)){
		MGFace* face2=new MGFace(*this);
		//std::cout<<(*face2)<<std::endl;///////////////*****:::::::::::
		for(int j=nboundaries-1; j>=0; j--){
			if(used_loops[j])
				face2->erase_boundary(j);
		}
		face2->add_boundary(new_loop.release());
		faces.push_back(face2);
	}
}

void set_used_loop(
	const std::deque<MGTrimLoop*>& used_tloops,
	std::vector<bool> used_loops
){
	size_t ntloop=used_tloops.size();
	for(size_t i=0; i<ntloop; i++){
		used_tloops[i]->set_used_loop_flag(used_loops);
	}
}

//Extract sub face that is bounded by network loops.
//Extracted sub face is the smallest closed part of this face bounded by
//the networks that includes the parameter position uv(u,v).
void MGFace::extract_sub_face(
	const MGPvector<MGLoop>& networks,
	const MGPosition& uv,
	std::auto_ptr<MGFace>& face//Result extracted face will be output.
)const{
	face=std::auto_ptr<MGFace>(new MGFace(*this));
	mgCPointsVec cpointsVec(*this,uv);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	MGPvector<MGLoop> iloops;
	int nnetworks=networks.size();
	for(int i=nnetworks-1; i>=0; i--){
		MGPvector<MGLoop> oloops;
		int closed_loop_obtained
			=cpointsVec.extract_trim_loop_info(*(networks[i]),oloops,iloops);
		if(closed_loop_obtained){
			face->add_boundary(oloops.release(0));
			return;
		}
	}

	std::auto_ptr<MGLoop> loop;
	std::deque<MGTrimLoop*> used_tloops;
	int eout=cpointsVec.extract_uv_boundary(loop,used_tloops);
	if(eout){
		size_t nboundaries=number_of_boundaries();
		std::vector<bool> used_loops(nboundaries);
		set_used_loop(used_tloops,used_loops);
		for(int j=nboundaries-1; j>=0; j--){
			if(used_loops[j])
				face->erase_boundary(j);
		}
		face->add_boundary(loop.release());
	}else{
		size_t niloops=iloops.size();
		for(size_t i=0; i<niloops; i++){
			face->add_boundary(iloops.release(i));
		}
		return;
	}
}

//Split the face giving networks loops. Splitting is done by finding the smallest
//closed areas out of networks.
void MGFace::split(
	const MGPvector<MGLoop>& networks,
	MGPvector<MGFace>& faces//Result trimmed face(s) will be appended.
)const{
	//std::cout<<networks<<std::endl;
	bool handled=false;
	mgCPointsVec cpointsVec(*this);
		//Found loops that are connected to
		//an inner loop of loop id loop_id will be stored
		//in this cpointsVec[loop_id].

	MGPvector<MGLoop> iloops,oloops;
	int nnetworks=networks.size();
	for(int ii=nnetworks-1; ii>=0; ii--)
		cpointsVec.extract_trim_loop_info(*(networks[ii]),oloops,iloops);
	//std::cout<<cpointsVec<<std::endl;

	size_t nolp=oloops.size();
	for(size_t i=0; i<nolp; i++){
		MGFace* nface=new MGFace(*this);
		nface->add_boundary(oloops.release(i));
		faces.push_back(nface);
		handled=true;
	}

	std::auto_ptr<MGFace> face2(new MGFace(*this));
	size_t nilp=iloops.size();
	for(size_t i=0; i<nilp; i++){
		face2->add_boundary(iloops.release(i));
		handled=true;
	}

	size_t nboundaries=number_of_boundaries();
	for(size_t i=0; i<nboundaries; i++){
		mgCPoints& lveci=cpointsVec[i];
		size_t nlv=lveci.size();
		assert((nlv%2)==0);//nlv must be even(pair of coming_in and going_out).
		if(!nlv)
			continue;

		for(size_t j=0; j<nlv; j++){
			std::auto_ptr<MGLoop> loop;
			std::deque<MGTrimLoop*> tloops;
			mgCONNECT_POINT& fpoint=lveci.m_cpoints[j];
			int extract_inf=cpointsVec.extract_loop(fpoint,loop,tloops);
			if(!extract_inf)
				continue;
		
			std::vector<bool> used_loops(nboundaries,false);
			size_t ntloop=tloops.size();
			for(size_t k=0; k<ntloop; k++){
				tloops[k]->set_used_loop_flag(used_loops);
				tloops[k]->set_null();
			}

			MGFace* nface2=face2->clone();
			for(int j=nboundaries-1; j>=0; j--){
				if(used_loops[j])
					nface2->erase_boundary(j);
			}
			nface2->add_boundary(loop.release());
			//std::cout<<(*nface2)<<std::endl;///////////////*****:::::::::::
			faces.push_back(nface2);
			handled=true;
		}
	}

	if(!handled)
		faces.push_back(face2.release());
}
