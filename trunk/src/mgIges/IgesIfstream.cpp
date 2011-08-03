/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Implementaion for class MGIgesIfstream.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Object.h"
#include "mg/CompositeCurve.h"
#include "mg/Cylinder.h"
#include "mg/Group.h"
#include "mgGL/LineStipple.h"
#include "mgGL/LineWidth.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD124.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mgGL/Color.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;

// Constructors.

// Creates an object of class MGIgesIfstream with a filename.
MGIgesIfstream::MGIgesIfstream(const char* file_name)
:m_ifstream(file_name){
	initialize();
}

void MGIgesIfstream::initialize(){
	if(!good())
		return;

	MGIgesFstream::initialize();

	char lineData[73];
	char ID;
	int sequence;

	//Construct m_StartSection, m_GSection, m_DirectryEntry, m_ParamDataLines.

	//1. Start Section
	get_one_line(lineData,ID,sequence);
	if(!good())
		return;

	while(ID=='S'){
		m_StartSection+=lineData;
		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}

	//2. Global Section
	std::string Gsecstring;
	m_nlineGSec=0;
	while(ID=='G'){
		Gsecstring+=lineData;
		m_nlineGSec++;
		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}
	m_GSection.read_in(Gsecstring);

	//3. Directory Entry Section
	MGIgesDirectoryEntry* de;
	while(ID=='D'){
		std::string DEstring(lineData);
		get_one_line(lineData,ID,sequence); assert(ID=='D');
		DEstring+=lineData;//Concatenate two DE lines into one.
		de=new MGIgesDirectoryEntry(DEstring);
		m_DirectoryEntries.push_back(de);

		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}

	char pDelimeter=m_GSection.paramDelimeter();//parameter delimeter

	//4. Parameter data Section
	while(ID=='P'){
		std::stringstream seq;seq.rdbuf()->str(lineData+65);
		//std::cout<<lineData<<std::endl;
		//std::cout<<seq.str()<<std::endl;
		int lnumber;seq>>lnumber;
		int DEpointer=MGIges::lnumber_to_DEpointer(lnumber);
		MGIgesDirectoryEntry& de=*(m_DirectoryEntries[DEpointer]);

		//Construct the string of the current Parameter data.
		int nlines=de.ParameterLineCount();
		std::string paramData(lineData,64);
		for(int i=1;i<nlines; i++){
			get_one_line(lineData,ID,sequence,64);
			paramData+=lineData;
			assert(ID=='P');
		}
		de.setPD(pDelimeter,paramData);

		//Get next line to process.
		get_one_line(lineData,ID,sequence);
		if(!good())
			return;
	}

	m_vertexListMap.set_ifstream(this);
	m_edgeListMap.set_ifstream(this);
}

//Open the file. When this is opened already, this will be closed, then will
//be opened.
void MGIgesIfstream::open(const char* file_name){
	if(is_open())
		close();
	m_ifstream.open(file_name);
	initialize();
}

void MGIgesIfstream::close(){
	clear();
	m_ifstream.close();
}

// Reads all objects of the IGES stored in the file
// as MGCL objects.
MGIgesIfstream& MGIgesIfstream::operator>>(MGGroup& group){
	int n=m_DirectoryEntries.size();//Num of directory entries.
	for(int i=1; i<n; i++){//i starts from 1 since 1st DE is dummy.
		MGIgesDirectoryEntry& de=*(m_DirectoryEntries[i]);
		if(!de.is_independent())
			continue;
		int tnum=de.EntityTypeNumber();
		if(tnum==TRANSFORMATION_MATRIX)//When transformation matrix
			continue;
		if(tnum==COLOR_DEFINITION)//When color definition
			continue;

		MGGel* gel=convert_to_gel(de);
		if(gel){
			//std::cout<<(*gel)<<std::endl;
			group.push_back(gel);
		}
	}
	return *this;
}

//From the current stream position, get one line data.
void MGIgesIfstream::get_one_line(
	char* lineData,	//line data without ID letter and sequence number
					//(that is, from column 1 to nchar) will be output.
					//buffer length must be >=(nchar+1).
	char& sectionID_letter,	//section identification letter of the line.
	int& sequence,	//ascending sequence number of the line.
	int nchar		//number of characters of one line
		//(When Parameter Data section nchar=64, and otherwise nchar=72)
){
	assert(nchar<=72);
	m_ifstream.get(lineData,nchar+1);
	char seqID[18];
	m_ifstream.get(seqID,74-nchar);
	sectionID_letter=seqID[72-nchar];
	m_ifstream>>sequence;//read sequence.
	char linefeed;
	m_ifstream.get(linefeed);	//read line feed.
}

//Convert i-th MGIgesDirectoryEntry object(m_DirectoryEntries[i])
//to MGObject that is a newed object.
//When de was not an independent object, null will be returned.
MGGel* MGIgesIfstream::convert_to_gel(
	int i
)const{
	MGIgesDirectoryEntry& de=*(m_DirectoryEntries[i]);
	return convert_to_gel(de);
}
MGGel* MGIgesIfstream::convert_to_gel(
	const MGIgesDirectoryEntry& de
)const{
	int typeNumber=de.EntityTypeNumber();
	MGGel* gel=0;
	switch(typeNumber){
	case CIRCULAR_ARC: gel=convert_arc(de); break;
	case COMPOSITE_CURVE: gel=convert_composite(de); break;
	case CONIC_ARC: gel=convert_conic_arc(de); break;
	case PLANE: gel=convert_plane(de); break;
	case LINE: gel=convert_line(de); break;
	case PARAMETRIC_SPLINE_CURVE: gel=convert_spline(de); break;
	case MGIges::POINT: gel=convert_point(de); break;
	case RULED_SURFACE: gel=convert_ruled_surface(de); break;
	case SURFACE_OF_REVOLUTION: gel=convert_revolution_surface(de); break;
	case TABULATED_CYLINDER: gel=convert_tab_cyl(de); break;
	case RATIONAL_BSPLINE_CURVE: gel=convert_nurbs(de);break;
	case RATIONAL_BSPLINE_SURFACE: gel=convert_nurbs_surface(de); break;
	case BOUNDED_SURFACE: gel=convert_bounded_surface(de); break;
	case TRIMMED_SURFACE: gel=convert_trimmed_surface(de); break;
	case SPHERE: gel=convert_sphere158(de); break;
	case MANIFOLD_SOLID_BREP_OBJECT: gel=convert_MSBO(de); break;
	case PLANE_SURFACE: gel=convert_planeSurface(de); break;
	case RIGHT_CIRCULAR_CYLINDRICAL_SURFACE: gel=convert_cylinder(de); break;
	case SPHERICAL_SURFACE: gel=convert_sphere(de); break;
	case ASSOCIATIVITY_INSTANCE: gel=convert_group(de); break;
	case FACE: gel=convert_face(de); break;
	case MGIges::SHELL: gel=convert_shell(de); break;
	default:std::cout<<"MGIgesIfstream::convert_to_gel:Non object typeNumber:"
				<<typeNumber<<std::endl;
	}

	if(!gel)
		return 0;

	transform(de,*gel);
	MGAttribedGel* agel=dynamic_cast<MGAttribedGel*>(gel);

	//Visibility.
	if(!de.is_visible()){
		agel->set_no_display();
	}

	//Color
	int color=de.ColorNumber();
	MGColor* mcolor=0;
	if(color>0){
		mcolor=new MGColor(MGColor::get_instance(static_cast<MGColor::ColorID>(color)));
	}else if(color<0){
		const MGIgesDirectoryEntry* color_de=directoryEntry(-color);
		mcolor=convert_color(*color_de);
	}
	if(mcolor)
		agel->set_GLattrib(mcolor);

	//Line width
	int lw=de.LineWeightNumber();
	if(lw)
		agel->set_GLattrib(new MGLineWidth(de.LineWidth(GSection())));

	//Line Font
	int lf=de.LineFontPattern();
	if(lf>1){
		agel->set_GLattrib(new MGLineStipple(MGLineStipple::LineFont(lf)));
	}

	//std::cout<<(*gel)<<std::endl;////***********::
	return gel;
}

//Transform obj if de has the transformation matrix.
void MGIgesIfstream::transform(
	const MGIgesDirectoryEntry& de,	//de of the object obj.
	MGGel& obj					//Object to transform.
)const{
	int tid=de.transformID();
	if(!tid)
		return;

	const MGIgesDirectoryEntry& trde=*(m_DirectoryEntries[tid]);
	const MGIgesPD124* pd124=static_cast<const MGIgesPD124*>(trde.paramData().get());
	MGTransf tr;
	pd124->convert_to_MGTransf(tr);
	obj.transform(tr);
}
