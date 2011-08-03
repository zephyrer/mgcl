/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD508.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/CompositeCurve.h"
#include "topo/BVertex.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD508.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD508 is the class for Iges parameter data type 508(LOOP Entity).

// Constructors.

//! Constructs an object of class MGIgesPD508.
MGIgesPD508::MGIgesPD508(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(LOOP,DEpointer), m_edges(1){
}

//Read in parameter data from string stream data.
void MGIgesPD508::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int num_edges;
	get_integer(pDelimeter,pdstream,num_edges);
	m_edges.resize(num_edges+1);
	for(int i=1; i<=num_edges; i++){
		MGIges508Edge* edge=new MGIges508Edge;
		edge->read_in(pDelimeter,pdstream);
		m_edges[i]=edge;
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD508::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	int num_edges=m_edges.size()-1;
	put_integer(num_edges,gsec,plines);
	for(int i=1; i<=num_edges; i++){
		const MGIges508Edge& edgei=*(m_edges[i]);
		edgei.write_out_into_string(gsec,plines);
	}
}

//Read in parameter data from string stream data.
void MGIges508Edge::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_type);
	get_DEpointer(pDelimeter,pdstream,m_edge_list);
	get_integer(pDelimeter,pdstream,m_edge);
	get_integer(pDelimeter,pdstream,m_orientation);
	int num_parameter_curve;
	get_integer(pDelimeter,pdstream,num_parameter_curve);
	m_isoparameterics.resize(num_parameter_curve+1);
	m_pcurves.resize(num_parameter_curve+1);
	for(int i=1; i<=num_parameter_curve; i++){
		int iso;
		get_integer(pDelimeter,pdstream,iso);
		m_isoparameterics[i]=iso ? 1:0;
		int pcurveDE;//pointer to the DE of the underlying parameter curve.
		get_DEpointer(pDelimeter,pdstream,pcurveDE);
		m_pcurves[i]=pcurveDE;
	}
}

//Write out this MGIges508Edge as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIges508Edge::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_integer(m_type,gsec,plines);
	put_DEpointer(m_edge_list,gsec,plines);
	put_integer(m_edge,gsec,plines);
	put_integer(m_orientation,gsec,plines);
	int num_parameter_curve=m_pcurves.size()-1;
	put_integer(num_parameter_curve,gsec,plines);
	for(int i=1; i<=num_parameter_curve; i++){
		int iso=m_isoparameterics[i] ? 1:0;
		put_integer(iso,gsec,plines);
		put_DEpointer(m_pcurves[i],gsec,plines);
	}
}

//Get the edge pointer(newed object).
MGEdge* MGIges508Edge::get_edge(const MGIgesIfstream& ifs)const{
	return ifs.m_edgeListMap.get_egde(m_edge_list,m_edge);
}

//Get the edge pointer(newed object).
MGBVertex* MGIges508Edge::get_BVertex(const MGIgesIfstream& ifs)const{
	return ifs.m_vertexListMap.get_BVertex(m_edge_list,m_edge);
}

//Convert de(type=508: LOOP entry) to MGLoop.
//Returned is a newed MGLoop object.
MGLoop* MGIgesIfstream::convert_loop(
	const MGIgesDirectoryEntry& de,	//directory entry of type 508.
	const MGSurface& srf //Base surface whose boundary this loop will make.
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD508* pd508=static_cast<const MGIgesPD508*>(pd.get());
	std::auto_ptr<MGLoop> loop=std::auto_ptr<MGLoop>(new MGLoop);
	size_t nedges=pd508->number_of_edges();
	if(!nedges)
		return 0;

	std::vector<double> pspan[2]; int peri_num[2];
	size_t startID=1;
	const MGIges508Edge& edge1=pd508->edge(1);
	MGEdge* bedge1;
	const MGCurve* wcurve1;
	if(!edge1.is_vertex()){
		bedge1=edge1.get_edge(*this);//This is binder edge.
		if(!bedge1)
			return 0;

		wcurve1=bedge1->base_curve();
		int nCom1st=srf.getPerimeterCommon(*wcurve1,pspan,peri_num);
		if(nCom1st>=2){
			startID=2;
		}else{
			double tLast=wcurve1->param_s(), terror=wcurve1->param_error();
			loop->append_edge_from_crv(srf,*wcurve1,tLast,terror,pspan[0],peri_num[0],
				edge1.direction_is_opposite());
			loop->last_edge()->set_binder(*bedge1);
		}
	}

	std::vector<double> pspan2[2]; int peri_num2[2];
	std::vector<double>& pspan0=pspan2[0];
	std::vector<double>& pspan1=pspan2[1];
	for(size_t i=2; i<=nedges; i++){
		const MGIges508Edge& edgei=pd508->edge(i);
		if(edgei.is_vertex())
			continue;

		MGEdge* bedge=edgei.get_edge(*this);//This is binder edge.
		if(!bedge)
			return 0;

		const MGCurve* wcurve=bedge->base_curve();
		int nComi=srf.getPerimeterCommon(*wcurve,pspan2,peri_num2);
		if(nComi>=2){
			MGPosition uv0=srf.perimeter_uv(peri_num2[0],pspan0[2]);
			MGPosition uv1=srf.perimeter_uv(peri_num2[1],pspan1[2]);
			MGPosition uvE=loop->end_point();
			if((uvE-uv1).len()<(uvE-uv0).len()){
				pspan0=pspan1; peri_num2[0]=peri_num2[1];
			}
		}
		double tLast=wcurve->param_s(), terror=wcurve->param_error();
		loop->append_edge_from_crv(srf,*wcurve,tLast,terror,pspan0,peri_num2[0],
			edgei.direction_is_opposite());
		loop->last_edge()->set_binder(*bedge);
	}
	if(startID==2){
		MGPosition uv0=srf.perimeter_uv(peri_num[0],pspan[0][2]);
		MGPosition uv1=srf.perimeter_uv(peri_num[1],pspan[1][2]);
		MGPosition uvE=loop->end_point();
		if((uvE-uv1).len()<(uvE-uv0).len()){
			pspan[0]=pspan[1]; peri_num[0]=peri_num[1];
		}
		double tLast=wcurve1->param_s(), terror=wcurve1->param_error();
		loop->append_edge_from_crv(srf,*wcurve1,tLast,terror,pspan[0],peri_num[0],
			edge1.direction_is_opposite());
		loop->last_edge()->set_binder(*bedge1);
	}
	loop->make_close();
	return loop.release();
}
