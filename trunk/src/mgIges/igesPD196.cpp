/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD196.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD116.h"
#include "mgiges/IgesPD123.h"
#include "mgiges/IgesPD196.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD196 is the class for Iges parameter data type 196(sphere surface).

// Constructors.

//! Constructs an object of class MGIgesPD196.
MGIgesPD196::MGIgesPD196(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(SPHERICAL_SURFACE,DEpointer),m_locationDE(0),m_radius(0.),m_axisDE(0),m_refdirDE(0){
}

//Construct PD196, supplying location on the plane, normal of the plane, and
//reference direction(might be omitted).
MGIgesPD196::MGIgesPD196(
	int locationDE, double radius, int axisDE, int refdirDE
):MGIgesPD(SPHERICAL_SURFACE),m_locationDE(locationDE),m_radius(radius),
m_axisDE(axisDE),m_refdirDE(refdirDE){
}

//Get the plane origin(LOCATION) into origin.
void MGIgesPD196::getCenter(const MGIgesIfstream& ifs, MGPosition& center)const{
	const MGIgesDirectoryEntry* de=ifs.directoryEntry(m_locationDE);
	assert(de->EntityTypeNumber()==116);
	if(de->EntityTypeNumber()!=116){
		center=MGPosition(0.,0.,0.);
		return;
	}
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD116* pd116=static_cast<const MGIgesPD116*>(pd.get());
	pd116->convert_to_position(center);
}

//Get the plane normal into nromal.
void MGIgesPD196::getAxis(const MGIgesIfstream& ifs, MGUnit_vector& axis)const{
	const MGIgesDirectoryEntry* de=ifs.directoryEntry(m_axisDE);
	assert(de->EntityTypeNumber()==DIRECTION);
	if(de->EntityTypeNumber()!=DIRECTION){
		axis=MGVector(0.,0.,1.);
		return;
	}
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD123* pd123=static_cast<const MGIgesPD123*>(pd.get());
	MGVector dir;
	pd123->convert_to_vector(dir);
	axis=dir;
}

//Get the plane reference direction(REFDIR) into refdir.
void MGIgesPD196::getRefdir(const MGIgesIfstream& ifs, MGVector& refdir)const{
	const MGIgesDirectoryEntry* de=ifs.directoryEntry(m_refdirDE);
	assert(de->EntityTypeNumber()==DIRECTION);
	if(de->EntityTypeNumber()!=DIRECTION){
		refdir=MGVector(1.,0.,0.);
		return;
	}
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD123* pd123=static_cast<const MGIgesPD123*>(pd.get());
	pd123->convert_to_vector(refdir);
}

//Read in parameter data from string stream data.
void MGIgesPD196::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_locationDE);
	get_real(pDelimeter,pdstream,m_radius);
	if(DEpointer()->FormNumber()==1){
		get_DEpointer(pDelimeter,pdstream,m_axisDE);
		get_DEpointer(pDelimeter,pdstream,m_refdirDE);
	}else{
		m_axisDE=0;
		m_refdirDE=0;
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD196::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_locationDE,gsec,plines);
	put_real(m_radius,gsec,plines);
	if(m_axisDE){
		put_DEpointer(m_axisDE,gsec,plines);
		put_DEpointer(m_refdirDE,gsec,plines);
	}
}

//Convert de to MGObject(a newed object). de must be of type SPHERICAL_SURFACE(196)
//(SPHERE surface).
//Output MGObject is a MGSphere.
MGSphere* MGIgesIfstream::convert_sphere(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	assert(de.EntityTypeNumber()==SPHERICAL_SURFACE);
	if(de.EntityTypeNumber()!=SPHERICAL_SURFACE)
		return 0;

	const MGIgesPD196* pd196=static_cast<const MGIgesPD196*>(pd.get());
	MGPosition center; pd196->getCenter(*this,center);
	double r=pd196->getRadius();
	if(de.FormNumber()==0)
		return new MGSphere(center,r);//nonparameterized sphere.

	//Parameterized.
	MGVector refdir;
	MGUnit_vector B;
	pd196->getAxis(*this,B);
	pd196->getRefdir(*this,refdir);
	return new MGSphere(center,r,B,refdir);//parameterized sphere.
}
