/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD192.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Cylinder.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD116.h"
#include "mgiges/IgesPD123.h"
#include "mgiges/IgesPD192.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD192 is the class for Iges parameter data type 192(cylinder surface).

// Constructors.

//! Constructs an object of class MGIgesPD192.
MGIgesPD192::MGIgesPD192(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(RIGHT_CIRCULAR_CYLINDRICAL_SURFACE,DEpointer),m_locationDE(0),m_normalDE(0),m_radius(0.),m_refdirDE(0){
}

//Construct PD192, supplying location on the plane, normal of the plane, and
//reference direction(might be omitted).
MGIgesPD192::MGIgesPD192(
	int locationDE, int normalDE, double radius, int refdirDE)
:MGIgesPD(RIGHT_CIRCULAR_CYLINDRICAL_SURFACE),m_locationDE(locationDE),m_normalDE(normalDE)
,m_radius(radius),m_refdirDE(refdirDE){
}

//Get the plane origin(LOCATION) into origin.
void MGIgesPD192::getOrigin(const MGIgesIfstream& ifs, MGPosition& origin)const{
	const MGIgesDirectoryEntry* de=ifs.directoryEntry(m_locationDE);
	assert(de->EntityTypeNumber()==116);
	if(de->EntityTypeNumber()!=116){
		origin=MGPosition(0.,0.,0.);
		return;
	}
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD116* pd116=static_cast<const MGIgesPD116*>(pd.get());
	pd116->convert_to_position(origin);
}

//Get the plane normal into nromal.
void MGIgesPD192::getNormal(const MGIgesIfstream& ifs, MGUnit_vector& normal)const{
	const MGIgesDirectoryEntry* de=ifs.directoryEntry(m_normalDE);
	assert(de->EntityTypeNumber()==DIRECTION);
	if(de->EntityTypeNumber()!=DIRECTION){
		normal=MGVector(0.,0.,1.);
		return;
	}
	const std::auto_ptr<MGIgesPD>& pd=de->paramData();
	const MGIgesPD123* pd123=static_cast<const MGIgesPD123*>(pd.get());
	MGVector dir;
	pd123->convert_to_vector(dir);
	normal=dir;
}

//Get the plane reference direction(REFDIR) into refdir.
void MGIgesPD192::getRefdir(const MGIgesIfstream& ifs, MGVector& refdir)const{
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
void MGIgesPD192::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_locationDE);
	get_DEpointer(pDelimeter,pdstream,m_normalDE);
	get_real(pDelimeter,pdstream,m_radius);
	if(DEpointer()->FormNumber()==1){
		get_DEpointer(pDelimeter,pdstream,m_refdirDE);
	}else
		m_refdirDE=0;
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD192::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_locationDE,gsec,plines);
	put_DEpointer(m_normalDE,gsec,plines);
	put_real(m_radius,gsec,plines);
	if(m_refdirDE)
		put_DEpointer(m_refdirDE,gsec,plines);
}

//Convert de to MGObject(a newed object). de must be of type RIGHT_CIRCULAR_CYLINDRICAL_SURFACE(192)
//(Right circular cylindrical surface).
//Output MGObject is a MGCylinder(unbounded infinite cylinder).
MGCylinder* MGIgesIfstream::convert_cylinder(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	assert(de.EntityTypeNumber()==192);
	if(de.EntityTypeNumber()!=192)
		return 0;

	const MGIgesPD192* pd192=static_cast<const MGIgesPD192*>(pd.get());
	MGPosition origin; pd192->getOrigin(*this,origin);
	MGUnit_vector normal; pd192->getNormal(*this,normal);

	MGVector refdir(1.,0.,0.);
	MGUnit_vector U, V;
	if(de.FormNumber()==1){//Parameterized.
		pd192->getRefdir(*this,refdir);
		U=refdir-(normal%refdir)*normal;
		V=normal*U;
	}else{
		normal.orthonormal(refdir,U,V);
	}
	double r=pd192->getRadius();
	MGInterval rng(0.,mgDBLPAI);
	MGEllipse circle(origin,U*r,V*r,rng);
	MGStraight axis(MGEReal(MGINFINITE_PLUS),MGEReal(MGINFINITE_MINUS),normal,origin);

	MGCylinder* cyl=new MGCylinder(circle,axis);
	return cyl;
}