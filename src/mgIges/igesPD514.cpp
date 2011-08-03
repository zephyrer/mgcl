/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD514.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD514.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD514 is the class for Iges parameter data type 514(SHELL).

// Constructors.

//! Constructs an object of class MGIgesPD514.
MGIgesPD514::MGIgesPD514(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(MGIges::SHELL,DEpointer),m_is_closed(false){
}

///append a face.
void MGIgesPD514::push_back(
		int face_DE,	///DE pointer of the face.
		bool same_direction
){
	m_faces.push_back(face_DE);
	m_orientations.push_back(same_direction);
}

//Read in parameter data from string stream data.
void MGIgesPD514::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	int num_faces;
	get_integer(pDelimeter,pdstream,num_faces);
	m_faces.resize(num_faces);
	m_orientations.resize(num_faces);
	for(int i=0; i<num_faces; i++){
		int faceDE;
		get_DEpointer(pDelimeter,pdstream,faceDE);
		m_faces[i]=faceDE;
		int orientation;
		get_integer(pDelimeter,pdstream,orientation);
		m_orientations[i]=orientation ? true:false;
	}
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD514::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	int num_faces=m_faces.size();
	put_integer(num_faces,gsec,plines);
	for(int i=0; i<num_faces; i++){
		put_DEpointer(m_faces[i],gsec,plines);
		int orientation=m_orientations[i] ? 1:0;
		put_integer(orientation,gsec,plines);
	}
}

//Convert de(type=514: SHELL) to MGShell.
//Returned is a newed MGShell object.
MGShell* MGIgesIfstream::convert_shell(
	const MGIgesDirectoryEntry& de
)const{
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD514* pd514=static_cast<const MGIgesPD514*>(pd.get());
	size_t nface=pd514->m_faces.size();
	if(!nface)
		return 0;

	MGShell* shell=0;
	for(size_t i=0; i<nface; i++){
		MGIgesDirectoryEntry& dei=*(m_DirectoryEntries[pd514->m_faces[i]]);
		MGFace* f=convert_face(dei);
		if(!f)
			continue;

		if(!pd514->m_orientations[i])
			f->negate();
		//std::cout<<*f<<std::endl;;///////*********
		if(shell)
			shell->append_face(f);
		else
			shell=new MGShell(f);
	}
	return shell;
}
