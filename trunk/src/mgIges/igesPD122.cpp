/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD122.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD122.h"
#include "mgiges/IgesGsec.h"
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD122 is the class for Iges parameter data type 122(Tabulated Cylinder).

// Constructors.

//! Constructs an object of class MGIgesPD122.
MGIgesPD122::MGIgesPD122(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(TABULATED_CYLINDER,DEpointer), m_directrix_DE(0){
	for(int i=0; i<3; i++)
		m_terminate_point[i]=0.;
}

//! Constructs an object of class MGIgesPD122.
MGIgesPD122::MGIgesPD122(int diretrix_DE, double terminate_point[3])
:MGIgesPD(TABULATED_CYLINDER), m_directrix_DE(diretrix_DE){
	for(int i=0; i<3; i++)
		m_terminate_point[i]=terminate_point[i];
}

//Read in parameter data from string stream data.
void MGIgesPD122::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_directrix_DE);
	get_real(pDelimeter,pdstream,m_terminate_point[0]);
	get_real(pDelimeter,pdstream,m_terminate_point[1]);
	get_real(pDelimeter,pdstream,m_terminate_point[2]);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD122::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_directrix_DE,gsec,plines);
	put_real(m_terminate_point[0],gsec,plines);
	put_real(m_terminate_point[1],gsec,plines);
	put_real(m_terminate_point[2],gsec,plines);

}

//Convert de(type=122: Ruled surface) to MGSurface.
//Returned is a newed object.
MGSurface* MGIgesIfstream::convert_tab_cyl(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD122* pd122=static_cast<const MGIgesPD122*>(pd.get());
	int directrixDE=pd122->m_directrix_DE;
	std::auto_ptr<MGGel> obj1=std::auto_ptr<MGGel>(convert_to_gel(directrixDE));
	
	std::auto_ptr<MGCurve> crv1;
	MGCurve* crv1p=dynamic_cast<MGCurve*>(obj1.get());
	if(crv1p)
		crv1=std::auto_ptr<MGCurve>(static_cast<MGCurve*>(obj1.release()));
	else
		return 0;

	MGVector generatrix=MGPosition(3,pd122->m_terminate_point) - crv1->start_point();
	std::auto_ptr<MGCurve> crv2(crv1->clone());
	(*crv2)+=generatrix;
	std::auto_ptr<MGSurface> surf=MGCL::create_ruled_surface(*crv1, *crv2);
	surf->change_range(true,0.,1.);
	surf->change_range(false,0.,1.);
	return surf.release();
}