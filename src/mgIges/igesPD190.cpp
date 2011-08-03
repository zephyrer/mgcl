/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD190.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesPD116.h"
#include "mgiges/IgesPD123.h"
#include "mgiges/IgesPD190.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD190 is the class for Iges parameter data type 100(plane surface).

// Constructors.

//! Constructs an object of class MGIgesPD190.
MGIgesPD190::MGIgesPD190(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(PLANE_SURFACE,DEpointer), m_locationDE(0), m_normalDE(0), m_refdirDE(0){
}

//Construct PD190, supplying location on the plane, normal of the plane, and
//reference direction(might be omitted).
MGIgesPD190::MGIgesPD190(
	int locationDE, int normalDE, int refdirDE
):MGIgesPD(PLANE_SURFACE), m_locationDE(locationDE), m_normalDE(normalDE), m_refdirDE(refdirDE){
}

//Get the plane origin(LOCATION) into origin.
void MGIgesPD190::getOrigin(const MGIgesIfstream& ifs, MGPosition& origin)const{
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
void MGIgesPD190::getNormal(const MGIgesIfstream& ifs, MGUnit_vector& normal)const{
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
void MGIgesPD190::getRefdir(const MGIgesIfstream& ifs, MGVector& refdir)const{
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
void MGIgesPD190::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	get_DEpointer(pDelimeter,pdstream,m_locationDE);
	get_DEpointer(pDelimeter,pdstream,m_normalDE);
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
void MGIgesPD190::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	put_DEpointer(m_locationDE,gsec,plines);
	put_DEpointer(m_normalDE,gsec,plines);
	if(m_refdirDE)
		put_DEpointer(m_refdirDE,gsec,plines);
}

//Convert de to MGObject(a newed object). de must be of type 190(Plane Surface).
//Output MGObject is a MGPlane(unbounded infinite plane).
MGPlane* MGIgesIfstream::convert_planeSurface(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	assert(de.EntityTypeNumber()==PLANE_SURFACE);
	if(de.EntityTypeNumber()!=PLANE_SURFACE)
		return 0;

	const MGIgesPD190* pd190=static_cast<const MGIgesPD190*>(pd.get());
	MGPosition origin; pd190->getOrigin(*this,origin);
	MGUnit_vector normal; pd190->getNormal(*this,normal);
	MGPlane* pl;
	if(de.FormNumber()==0)//Unparameterized.
		pl=new MGPlane(normal,origin);
	else{//Parameterized.
		MGVector refdir;
		pd190->getRefdir(*this,refdir);
		MGVector U=refdir-(normal%refdir)*normal;
		MGVector V=normal*U;
		pl=new MGPlane(U,V,origin);
	}
	return pl;
}

