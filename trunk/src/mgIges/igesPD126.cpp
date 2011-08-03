/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD126.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Plane.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD126.h"
#include "mgiges/IgesGsec.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD126 is the class for Iges parameter data type 126(NURBS).

// Constructors.

//! Constructs an object of class MGIgesPD126.
MGIgesPD126::MGIgesPD126(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(RATIONAL_BSPLINE_CURVE,DEpointer),m_upper_index(-1),m_degree(0),
m_planar(0),m_closed(0),m_non_rational(0),m_periodic(0),
m_start_param(0.),m_end_param(0.){
	m_normal[0]=0.;m_normal[1]=0.;m_normal[2]=1.;
}

//! Constructs an object of class MGIgesPD126.
MGIgesPD126::MGIgesPD126(const MGLBRep& lb)
:MGIgesPD(RATIONAL_BSPLINE_CURVE),m_upper_index(lb.bdim()-1),m_degree(lb.order()-1),
m_planar(0),m_closed(0),m_non_rational(1),m_periodic(0),
m_knots(lb.knot_vector()),m_weights(lb.bdim(),1.),
m_control_points(lb.line_bcoef()),
m_start_param(lb.param_s()),m_end_param(lb.param_e()){
	MGPlane pl; MGStraight sl; MGPosition P;
	if(lb.planar(pl,sl,P))
		m_planar=1;
	if(lb.is_closed())
		m_closed=1;
	const MGUnit_vector& N=pl.normal();
	m_normal[0]=N[0];m_normal[1]=N[1];m_normal[2]=N[2];
}

//! Constructs an object of class MGIgesPD126.
MGIgesPD126::MGIgesPD126(const MGRLBRep& lb)
:MGIgesPD(RATIONAL_BSPLINE_CURVE),m_upper_index(lb.bdim()-1),m_degree(lb.order()-1),
m_planar(0),m_closed(0),m_non_rational(0),m_periodic(0),
m_knots(lb.knot_vector()),m_weights(lb.bdim(),1.),
m_control_points(lb.bdim(),lb.sdim()),
m_start_param(lb.param_s()),m_end_param(lb.param_e()){
	MGPlane pl; MGStraight sl; MGPosition P;
	if(lb.planar(pl,sl,P))
		m_planar=1;
	if(lb.is_closed())
		m_closed=1;
	const MGUnit_vector& N=pl.normal();
	m_normal[0]=N[0];m_normal[1]=N[1];m_normal[2]=N[2];

	int n=lb.bdim();
	const MGBPointSeq& c=lb.line_bcoef();
	int sd=lb.sdim();
	for(int i=0; i<n; i++){
		double wi=c(i,sd);
		m_weights[i]=wi;
		for(int j=0; j<sd; j++)
			m_control_points(i,j)=c(i,j)/wi;
	}
}

//Read in parameter data from string stream data.
void MGIgesPD126::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_upper_index);
	get_integer(pDelimeter,pdstream,m_degree);
	get_integer(pDelimeter,pdstream,m_planar);
	get_integer(pDelimeter,pdstream,m_closed);
	get_integer(pDelimeter,pdstream,m_non_rational);
	get_integer(pDelimeter,pdstream,m_periodic);

	int order=m_degree+1;
	int nBrep=m_upper_index+1;

	//Read knot vector.
	int nBrepPorder=nBrep+order, i;
	m_knots.size_change(order,nBrep);
	for(i=0; i<nBrepPorder; i++)
		get_real(pDelimeter,pdstream,m_knots[i]);
	//if(m_knots[m_degree]>=m_knots[nBrep])
	//	std::cout<<m_knots<<std::endl;

	//Read weights.
	m_weights.resize(nBrep);
	for(i=0; i<nBrep; i++)
		get_real(pDelimeter,pdstream,m_weights[i]);

	//Read control points.
	m_control_points.resize(nBrep,3);
	for(i=0; i<nBrep; i++){
		get_real(pDelimeter,pdstream,m_control_points(i,0));
		get_real(pDelimeter,pdstream,m_control_points(i,1));
		get_real(pDelimeter,pdstream,m_control_points(i,2));
	}

	get_real(pDelimeter,pdstream,m_start_param);
	get_real(pDelimeter,pdstream,m_end_param);

	//Read normal.
	get_real(pDelimeter,pdstream,m_normal[0]);
	get_real(pDelimeter,pdstream,m_normal[1]);
	get_real(pDelimeter,pdstream,m_normal[2]);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD126::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_integer(m_upper_index,gsec,plines);
	put_integer(m_degree,gsec,plines);
	put_integer(m_planar,gsec,plines);
	put_integer(m_closed,gsec,plines);
	put_integer(m_non_rational,gsec,plines);
	put_integer(m_periodic,gsec,plines);

	int i;

	//Write out knot vector.
	int n=m_knots.length();
	for(i=0; i<n; i++)
		put_real(m_knots[i],gsec,plines);

	//Write out weights.
	int m=m_weights.size();
	for(i=0; i<m; i++)
		put_real(m_weights[i],gsec,plines);

	//Write out control points.
	int l=m_control_points.length();
	for(i=0; i<l; i++){
		put_real(m_control_points(i,0),gsec,plines);
		put_real(m_control_points(i,1),gsec,plines);
		put_real(m_control_points(i,2),gsec,plines);
	}
	put_real(m_start_param,gsec,plines);
	put_real(m_end_param,gsec,plines);

	//Write out normal.
	put_real(m_normal[0],gsec,plines);
	put_real(m_normal[1],gsec,plines);
	put_real(m_normal[2],gsec,plines);
}

//Convert de(type=126: Rational b-spline) to MGCurve(MGLBRep or MBRLBRep).
//Returned is a newed object.
MGCurve* MGIgesIfstream::convert_nurbs(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD126* pd126=static_cast<const MGIgesPD126*>(pd.get());
	MGCurve* crv;
	if(pd126->m_non_rational){
		crv=new MGLBRep(pd126->m_knots,pd126->m_control_points);
	}else{
		crv=new MGRLBRep(pd126->m_knots,pd126->m_control_points,pd126->m_weights);
	}
	//std::cout<<(*crv)<<std::endl;///////************
	return crv;
}
