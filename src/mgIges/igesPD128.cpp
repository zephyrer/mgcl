/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD128.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD128.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesPD128 is the class for Iges parameter data type 128(NURBS Surface).
using namespace MGIges;

// Constructors.

//! Constructs an object of class MGIgesPD128.
MGIgesPD128::MGIgesPD128(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(RATIONAL_BSPLINE_SURFACE,DEpointer),m_upper_indexU(-1),m_upper_indexV(-1),
m_degreeU(0),m_degreeV(0),m_closedU(0),m_closedV(0),m_non_rational(0),
m_periodicU(0),m_periodicV(0),m_start_paramU(0.),m_end_paramU(0.),
m_start_paramV(0.),m_end_paramV(0.){
}

//! Constructs an object of class MGIgesPD128.
MGIgesPD128::MGIgesPD128(const MGSBRep& sb)
:MGIgesPD(RATIONAL_BSPLINE_SURFACE),m_upper_indexU(sb.bdim_u()-1),m_upper_indexV(sb.bdim_v()-1),
m_degreeU(sb.order_u()-1),m_degreeV(sb.order_v()-1),m_closedU(0),m_closedV(0),
m_non_rational(1),m_periodicU(0),m_periodicV(0),
m_knotsU(sb.knot_vector_u()),m_knotsV(sb.knot_vector_v()),
m_weights(sb.bdim_u(),sb.bdim_v(),1),
m_control_points(sb.surface_bcoef()),
m_start_paramU(sb.param_s_u()),m_end_paramU(sb.param_e_u()),
m_start_paramV(sb.param_s_v()),m_end_paramV(sb.param_e_v()){
	int nu=sb.bdim_u(), nv=sb.bdim_v();
	for(int i=0; i<nu; i++){
		for(int j=0; j<nv; j++){
			m_weights(i,j,0)=1.;
		}
	}
}

//! Constructs an object of class MGIgesPD128(NURBS Surface)..
MGIgesPD128::MGIgesPD128(const MGRSBRep& sb)
:MGIgesPD(RATIONAL_BSPLINE_SURFACE),m_upper_indexU(sb.bdim_u()-1),m_upper_indexV(sb.bdim_v()-1),
m_degreeU(sb.order_u()-1),m_degreeV(sb.order_v()-1),m_closedU(0),m_closedV(0),
m_non_rational(0),m_periodicU(0),m_periodicV(0),
m_knotsU(sb.knot_vector_u()),m_knotsV(sb.knot_vector_v()),
m_weights(sb.bdim_u(),sb.bdim_v(),1),
m_control_points(sb.bdim_u(),sb.bdim_v(),sb.sdim()),
m_start_paramU(sb.param_s_u()),m_end_paramU(sb.param_e_u()),
m_start_paramV(sb.param_s_v()),m_end_paramV(sb.param_e_v()){
	int nu=sb.bdim_u(), nv=sb.bdim_v();
	const MGSPointSeq& c=sb.surface_bcoef();
	int sd=sb.sdim();
	for(int i=0; i<nu; i++){
		for(int j=0; j<nv; j++){
			double wij=c(i,j,sd);
			m_weights(i,j,0)=wij;
			for(int k=0; k<sd; k++)
				m_control_points(i,j,k)=c(i,j,k)/wij;
		}
	}
}

//Read in parameter data from string stream data.
void MGIgesPD128::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_upper_indexU);
	get_integer(pDelimeter,pdstream,m_upper_indexV);
	get_integer(pDelimeter,pdstream,m_degreeU);
	get_integer(pDelimeter,pdstream,m_degreeV);

	get_integer(pDelimeter,pdstream,m_closedU);
	get_integer(pDelimeter,pdstream,m_closedV);
	get_integer(pDelimeter,pdstream,m_non_rational);
	get_integer(pDelimeter,pdstream,m_periodicU);
	get_integer(pDelimeter,pdstream,m_periodicV);

	int orderU=m_degreeU+1, orderV=m_degreeV+1;
	int nBrepU=m_upper_indexU+1, nBrepV=m_upper_indexV+1;

	//Read knot vector.
	int nBrepUPorderU=nBrepU+orderU, i;
	m_knotsU.size_change(orderU,nBrepU);
	for(i=0; i<nBrepUPorderU; i++)
		get_real(pDelimeter,pdstream,m_knotsU[i]);

	int nBrepVPorderV=nBrepV+orderV, j;
	m_knotsV.size_change(orderV,nBrepV);
	for(j=0; j<nBrepVPorderV; j++)
		get_real(pDelimeter,pdstream,m_knotsV[j]);
	
	int nBrepUBYnBrepV=nBrepU*nBrepV;

	//Read weights.
	m_weights.resize(nBrepU,nBrepV,1);
	for(j=0; j<nBrepV; j++){
		for(i=0; i<nBrepU; i++)
			get_real(pDelimeter,pdstream,m_weights(i,j,0));
	}

	//Read control points..
	m_control_points.resize(nBrepU,nBrepV,3);
	for(j=0; j<nBrepV; j++){
		for(i=0; i<nBrepU; i++){
			get_real(pDelimeter,pdstream,m_control_points(i,j,0));
			get_real(pDelimeter,pdstream,m_control_points(i,j,1));
			get_real(pDelimeter,pdstream,m_control_points(i,j,2));
		}
	}

	get_real(pDelimeter,pdstream,m_start_paramU);
	get_real(pDelimeter,pdstream,m_end_paramU);
	get_real(pDelimeter,pdstream,m_start_paramV);
	get_real(pDelimeter,pdstream,m_end_paramV);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD128::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_integer(m_upper_indexU,gsec,plines);
	put_integer(m_upper_indexV,gsec,plines);
	put_integer(m_degreeU,gsec,plines);
	put_integer(m_degreeV,gsec,plines);
	put_integer(m_closedU,gsec,plines);
	put_integer(m_closedV,gsec,plines);
	put_integer(m_non_rational,gsec,plines);
	put_integer(m_periodicU,gsec,plines);
	put_integer(m_periodicV,gsec,plines);
	int knotsUlength=m_knotsU.length();

	//Write out knot vector.
	int i,j,k;
	for(i=0;i<knotsUlength;i++){
		put_real(m_knotsU[i],gsec,plines);
	}
	int knotsVlength=m_knotsV.length();
	for(i=0;i<knotsVlength;i++){
		put_real(m_knotsV[i],gsec,plines);
	}
	int nBrepU=m_upper_indexU+1, nBrepV=m_upper_indexV+1;

	//Write out weights.
	for(j=0;j<nBrepV;j++){
		for(i=0;i<nBrepU;i++){
			put_real(m_weights(i,j,0),gsec,plines);
		}
	}

	//Write out control points..
	for(j=0;j<nBrepV;j++){
		for(i=0;i<nBrepU;i++){
			for(k=0;k<3;k++){
				put_real(m_control_points(i,j,k),gsec,plines);
			}
		}
	}

	put_real(m_start_paramU,gsec,plines);
	put_real(m_end_paramU,gsec,plines);
	put_real(m_start_paramV,gsec,plines);
	put_real(m_end_paramV,gsec,plines);
}

//Convert de(type=128: Rational b-spline surface) to MGSurface(MGSBRep or MGRSBRep).
//Returned is a newed object.
MGSurface* MGIgesIfstream::convert_nurbs_surface(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD128* pd128=static_cast<const MGIgesPD128*>(pd.get());
	MGSurface* srf;
	if(pd128->m_non_rational){
		MGSBRep* sb=new MGSBRep(pd128->m_control_points,pd128->m_knotsU,pd128->m_knotsV);
		srf=sb;
//		if(sb->bdim_u()==7 && sb->bdim_v()==43){/////////////**********DEBUGWRITE
//			MGSPointSeq& sp=sb->surface_bcoef();
//			MGKnotVector& tu=sb->knot_vector_u();
//			MGKnotVector& tv=sb->knot_vector_v();
//			if(1.840<tu(7) && tu(0)<1.841 && 4.705<tv(43) && tv(43)<4.706
//				&& 165.296<sp(0,0,0)&& 165.298>sp(0,0,0))
//				std::cout<<(*sb)<<std::endl;
//		}/////////////**********DEBUGWRITE
	}else{
		srf=new MGRSBRep(pd128->m_control_points,pd128->m_weights,
			pd128->m_knotsU,pd128->m_knotsV);
	}
	//std::cout<<*srf<<std::endl;;///////*********
	return srf;
}
