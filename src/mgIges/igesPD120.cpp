/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD120.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD120.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD120 is the class for Iges parameter data type 120(Surface of Revolution).

// Constructors.

//! Constructs an object of class MGIgesPD120.
MGIgesPD120::MGIgesPD120(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(SURFACE_OF_REVOLUTION,DEpointer), m_axis_of_revolution_DE(0),m_generatrix_DE(0),
m_start_angle(0.), m_terminate_angle(0.){
}

//Read in parameter data from string stream data.
void MGIgesPD120::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_axis_of_revolution_DE);
	get_DEpointer(pDelimeter,pdstream,m_generatrix_DE);
	get_real(pDelimeter,pdstream,m_start_angle);
	get_real(pDelimeter,pdstream,m_terminate_angle);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD120::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_axis_of_revolution_DE,gsec,plines);
	put_DEpointer(m_generatrix_DE,gsec,plines);
	put_real(m_start_angle,gsec,plines);
	put_real(m_terminate_angle,gsec,plines);
}

//Convert de(type=120: Ruled surface) to MGSurface.
//Returned is a newed object.
MGSurface* MGIgesIfstream::convert_revolution_surface(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD120* pd120=static_cast<const MGIgesPD120*>(pd.get());
	int axisDE=pd120->m_axis_of_revolution_DE;
	std::auto_ptr<MGGel> obj1=std::auto_ptr<MGGel>(convert_to_gel(axisDE));
	MGStraight* axisp=dynamic_cast<MGStraight*>(obj1.get());
	if(!axisp)
		return 0;
	std::auto_ptr<MGStraight> axis=std::auto_ptr<MGStraight>(static_cast<MGStraight*>(obj1.release()));

	int generatrixDE=pd120->m_generatrix_DE;
	std::auto_ptr<MGGel> obj2=std::auto_ptr<MGGel>(convert_to_gel(generatrixDE));
	MGCurve* generatrixp=dynamic_cast<MGCurve*>(obj2.get());
	if(!generatrixp)
		return 0;

	std::auto_ptr<MGCurve> generatrix=std::auto_ptr<MGCurve>(static_cast<MGCurve*>(obj2.release()));
	MGTransf tr;
	tr.set_rotate_3D(axis->direction(),pd120->m_start_angle,axis->root_point());
	(*generatrix)*=tr;

	double sa=pd120->m_start_angle, ta=pd120->m_terminate_angle;
	double angle=ta - sa;
	std::auto_ptr<MGSurface> surf=MGCL::create_revolved_surface(*generatrix, *axis, angle);
	surf->exchange_uv();surf->negate(0);
	surf->change_range(false,sa,ta);
	return surf.release();
}
