/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Implementation for class MGIgesPD112(Spline curve).
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/PPRep.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD112.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD112 is the class for Iges parameter data type 112(LIne).

// Constructors.

//! Constructs an object of class MGIgesPD112.
MGIgesPD112::MGIgesPD112(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(PARAMETRIC_SPLINE_CURVE,DEpointer){
}

//Read in parameter data from string stream data.
void MGIgesPD112::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_integer(pDelimeter,pdstream,m_spline_type);
	get_integer(pDelimeter,pdstream,m_continuity);
	get_integer(pDelimeter,pdstream,m_dimension);

	int number_of_segments;
	get_integer(pDelimeter,pdstream,number_of_segments);
	int j;
	int number_of_bp=number_of_segments+1;//Number of break points.
	m_tau.resize(number_of_bp);
	for(j=0; j<number_of_bp; j++)
		get_real(pDelimeter,pdstream,m_tau[j]);//Read in break points.

	m_coefs.resize(number_of_bp);
	for(j=0; j<number_of_bp; j++){//j for break points.
		MGIgesSplineCoef& abcd=m_coefs[j];
		for(int k=0; k<3; k++){//k for (x,y,z)=space dimension
			for(int i=0; i<4; j++)//i for (A,B,C,D)=order
				get_real(pDelimeter,pdstream,abcd.m_abcd[i][k]);//Read in coefficients.
		}
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than 2.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD112::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_integer(m_spline_type,gsec,plines);
	put_integer(m_continuity,gsec,plines);
	put_integer(m_dimension,gsec,plines);

	int number_of_bp=m_tau.length();
	int number_of_segments=number_of_bp-1;
	put_integer(number_of_segments,gsec,plines);
	int i;
	for(i=0; i<number_of_bp; i++){
		put_real(m_tau[i],gsec,plines);
	}
	int csize=m_coefs.size();
	for(i=0; i<csize; i++){
		MGIgesSplineCoef abcd=m_coefs[i];//i番目のm_coefsクラスをインスタンス化
		for(int j=0; j<3; j++){
			for(int k=0; k<4; k++){
				put_real(abcd.m_abcd[k][j],gsec,plines);
			}
		}
	}
}

//Convert de to MGObject(a newed object). de must be of type 112(conventional cubic spline).
//Output MGObject is an MGLBRep object.
MGLBRep* MGIgesIfstream::convert_spline(
	const MGIgesDirectoryEntry& de
)const{
	int fnum=de.FormNumber();
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD112* pd112=static_cast<const MGIgesPD112*>(pd.get());
	int order=4;
	if(pd112->m_spline_type ==1)//linear
		order=2;
	else if(pd112->m_spline_type ==2)//quardric
		order=3;

	const MGNDDArray& tau=pd112->m_tau;//Break point seq.
	MGPPRep pprep(order, 3, tau);//Space dimension is 3.
	int num_seg=tau.length()-1;

	const double factorial[4]={1.,1.,2.,6.};
	for(int j=0; j<num_seg; j++){//j for break points.
		const MGIgesSplineCoef& abcd=pd112->m_coefs[j];
		for(int k=0; k<3; k++){//k for (x,y,z)=space dimension
			for(int i=0; i<order; j++){//i for (A,B,C,D)=order
				pprep.coef(i,j,k)=factorial[i]*abcd.m_abcd[i][k];
			}
		}
	}

	MGLBRep* lb=new MGLBRep(pprep);
	return lb;
}
