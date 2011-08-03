/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesDirectoryEntry.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Point.h"
#include "mg/Group.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/CompositeCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "topo/Face.h"
#include "mgiges/IgesDirectoryEntry.h"
#include "mgiges/IgesPD100.h"
#include "mgiges/IgesPD102.h"
#include "mgiges/IgesPD104.h"
#include "mgiges/IgesPD108.h"
#include "mgiges/IgesPD110.h"
#include "mgiges/IgesPD112.h"
#include "mgiges/IgesPD116.h"
#include "mgiges/IgesPD118.h"
#include "mgiges/IgesPD120.h"
#include "mgiges/IgesPD122.h"
#include "mgiges/IgesPD123.h"
#include "mgiges/IgesPD124.h"
#include "mgiges/IgesPD126.h"
#include "mgiges/IgesPD128.h"
#include "mgiges/IgesPD186.h"
#include "mgiges/IgesPD141.h"
#include "mgiges/IgesPD142.h"
#include "mgiges/IgesPD143.h"
#include "mgiges/IgesPD144.h"
#include "mgiges/IgesPD158.h"
#include "mgiges/IgesPD190.h"
#include "mgiges/IgesPD192.h"
#include "mgiges/IgesPD196.h"
#include "mgiges/IgesPD314.h"
#include "mgiges/IgesPD402.h"
#include "mgiges/IgesPD502.h"
#include "mgiges/IgesPD504.h"
#include "mgiges/IgesPD508.h"
#include "mgiges/IgesPD510.h"
#include "mgiges/IgesPD514.h"

using namespace std;
using namespace MGIges;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

MGEXTERN const MGIgesDEStatusNumber mgIgesDEStatusNumber=MGIgesDEStatusNumber();

//!	@brief MGIgesDirectoryEntry describes a directory entry section of an IGES file.
//In IGES file one directory entry is stored in two records. MGIgesDirectoryEntry
//holds the two records data.

//! Constructs an object of class MGIgesDirectoryEntry.
//Default constructor, includes all the defalut value of MGCL.
MGIgesDirectoryEntry::MGIgesDirectoryEntry()
:m_EntityTypeNumber(0), m_ParameterDataLine(0), m_Structure(0), m_LineFontPattern(0)
, m_Level(0), m_View(0), m_TransformationMatrix(0), m_LabelDisplayAssociativity(0)
, m_LineWeightNumber(0), m_ColorNumber(0), m_ParameterLineCount(0)
, m_FormNumber(0), m_EntitySubscriptNumber(0), m_ParamData(0){;}

//Construct from a string that contais one pair of lines of IGES File.
MGIgesDirectoryEntry::MGIgesDirectoryEntry(
	const string& DEstring
):m_EntityTypeNumber(0), m_ParameterDataLine(0), m_Structure(0), m_LineFontPattern(0)
, m_Level(0), m_View(0), m_TransformationMatrix(0), m_LabelDisplayAssociativity(0)
, m_LineWeightNumber(0), m_ColorNumber(0), m_ParameterLineCount(0)
, m_FormNumber(0), m_EntitySubscriptNumber(0), m_ParamData(0){
	stringstream(DEstring.substr(0,8))>>skipws>>m_EntityTypeNumber;
	stringstream(DEstring.substr(8,8))>>skipws>>m_ParameterDataLine;
	stringstream(DEstring.substr(16,8))>>skipws>>m_Structure;
	m_Structure=MGIges::lnumber_to_DEpointer(m_Structure);

	stringstream(DEstring.substr(24,8))>>skipws>>m_LineFontPattern;
	stringstream(DEstring.substr(32,8))>>skipws>>m_Level;
	stringstream(DEstring.substr(40,8))>>skipws>>m_View;
	m_View=MGIges::lnumber_to_DEpointer(m_View);

	stringstream(DEstring.substr(48,8))>>skipws>>m_TransformationMatrix;
	m_TransformationMatrix=MGIges::lnumber_to_DEpointer(m_TransformationMatrix);

	stringstream(DEstring.substr(56,8))>>skipws>>m_LabelDisplayAssociativity;
	m_LabelDisplayAssociativity=MGIges::lnumber_to_DEpointer(m_LabelDisplayAssociativity);
	m_StatusNumber.read_in(DEstring.substr(64,8));
	stringstream(DEstring.substr(80,8))>>skipws>>m_LineWeightNumber;
	stringstream(DEstring.substr(88,8))>>skipws>>m_ColorNumber;
	if(m_ColorNumber<0)
		m_ColorNumber=-MGIges::lnumber_to_DEpointer(-m_ColorNumber);

	stringstream(DEstring.substr(96,8))>>skipws>>m_ParameterLineCount;
	stringstream(DEstring.substr(104,8))>>skipws>>m_FormNumber;
	//stringstream(DEstring.substr(112,8))>>skipws>>m_Reserved1;
	//stringstream(DEstring.substr(120,8))>>skipws>>m_Reserved2;
	m_EntityLabel=DEstring.substr(128,8);
	stringstream(DEstring.substr(136,8))>>skipws>>m_EntitySubscriptNumber;
}
MGIgesDirectoryEntry::MGIgesDirectoryEntry(
	int EntityTypeNumber,
	const string& EntityLabel,
	const MGIgesDEStatusNumber& StatusNumber,
	int ColorNumber,//When negated value, is a pointer to the directory entry.
	int LineWeightNumber,
	int LineFontPattern,//When negated value, is a pointer to the directory entry.
	int TRANSFORMATION_MATRIX,//Pointer to the directory entity of the matrix.
	int FormNumber,
	int EntitySubscriptNumber
):m_EntityTypeNumber(EntityTypeNumber), m_ParameterDataLine(0), m_Structure(0),
m_LineFontPattern(LineFontPattern), m_Level(0), m_View(0),
m_TransformationMatrix(TRANSFORMATION_MATRIX), m_LabelDisplayAssociativity(0),
m_StatusNumber(StatusNumber), m_LineWeightNumber(LineWeightNumber),
m_ColorNumber(ColorNumber), m_ParameterLineCount(0), m_FormNumber(FormNumber),
m_EntityLabel(EntityLabel), m_EntitySubscriptNumber(EntitySubscriptNumber),
m_ParamData(0){;}

//Get line width.
float MGIgesDirectoryEntry::LineWidth(const MGIgesGSec& gsec)const{
	float w=float(gsec.m_width_of_max_line_weight);
	float mnum=float(gsec.m_max_number_of_line_weight_gradations);
	float lwnum=float(LineWeightNumber());
	return lwnum*w/mnum;
}

//Set pd in m_ParamData.
void MGIgesDirectoryEntry::setParamData(auto_ptr<MGIgesPD>& pd){
	assert(pd->DEpointer()==0);
	m_ParamData=pd;
	m_ParamData->m_DEpointer=this;
}

void MGIgesDirectoryEntry::setPD(auto_ptr<MGIgesPD>& pd){
	pd->setDE(this);
	m_ParamData=pd;
}

//Convert paramData to the each type parameter data, and set it in m_ParamData.
void MGIgesDirectoryEntry::setPD(
	char pDelimeter,		//parameter delimeter
	const string& paramData	//whole parameter data string.
							//Each of the parameter is separated by pDelimeter.
){
	istringstream pdstream(paramData);
	int typeNumber;
	pdstream>>typeNumber;
	string dummy;
	getline(pdstream,dummy,pDelimeter);
	assert(typeNumber==m_EntityTypeNumber);

	auto_ptr<MGIgesPD> pd;
	switch(typeNumber){
		case CIRCULAR_ARC: pd=auto_ptr<MGIgesPD>(new MGIgesPD100); break;
		case COMPOSITE_CURVE: pd=auto_ptr<MGIgesPD>(new MGIgesPD102); break;
		case CONIC_ARC: pd=auto_ptr<MGIgesPD>(new MGIgesPD104); break;
		case PLANE: pd=auto_ptr<MGIgesPD>(new MGIgesPD108); break;
		case LINE: pd=auto_ptr<MGIgesPD>(new MGIgesPD110); break;
		case PARAMETRIC_SPLINE_CURVE: pd=auto_ptr<MGIgesPD>(new MGIgesPD112); break;
		case MGIges::POINT: pd=auto_ptr<MGIgesPD>(new MGIgesPD116); break;
		case RULED_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD118); break;
		case SURFACE_OF_REVOLUTION: pd=auto_ptr<MGIgesPD>(new MGIgesPD120); break;
		case TABULATED_CYLINDER: pd=auto_ptr<MGIgesPD>(new MGIgesPD122); break;
		case DIRECTION: pd=auto_ptr<MGIgesPD>(new MGIgesPD123); break;
		case TRANSFORMATION_MATRIX: pd=auto_ptr<MGIgesPD>(new MGIgesPD124); break;
		case RATIONAL_BSPLINE_CURVE: pd=auto_ptr<MGIgesPD>(new MGIgesPD126); break;
		case RATIONAL_BSPLINE_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD128); break;
		case BOUNDARY: pd=auto_ptr<MGIgesPD>(new MGIgesPD141); break;
		case CURVE_ON_PARAMETRIC_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD142); break;
		case BOUNDED_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD143); break;
		case TRIMMED_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD144); break;
		case SPHERE: pd=auto_ptr<MGIgesPD>(new MGIgesPD158); break;
		case MANIFOLD_SOLID_BREP_OBJECT: pd=auto_ptr<MGIgesPD>(new MGIgesPD186); break;
		case PLANE_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD190); break;
		case RIGHT_CIRCULAR_CYLINDRICAL_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD192); break;
		case SPHERICAL_SURFACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD196); break;
		case COLOR_DEFINITION: pd=auto_ptr<MGIgesPD>(new MGIgesPD314); break;
		case ASSOCIATIVITY_INSTANCE: pd=auto_ptr<MGIgesPD>(new MGIgesPD402); break;
		case VERTEX: pd=auto_ptr<MGIgesPD>(new MGIgesPD502); break;
		case EDGE: pd=auto_ptr<MGIgesPD>(new MGIgesPD504); break;
		case LOOP: pd=auto_ptr<MGIgesPD>(new MGIgesPD508); break;
		case FACE: pd=auto_ptr<MGIgesPD>(new MGIgesPD510); break;
		case MGIges::SHELL: pd=auto_ptr<MGIgesPD>(new MGIgesPD514); break;
		default:
			std::cout<<"MGIgesDirectoryEntry::setPD:Usupported typeNumber:"
				<<typeNumber<<std::endl;
	}

	MGIgesPD* pdptr=pd.get();
	if(pdptr){
		setParamData(pd);
		pdptr->read_in(pDelimeter,pdstream);
	}
}

//Put this DE data to 2 lines of IGES file string data.
void MGIgesDirectoryEntry::put_to_string(
	int ParameterDataLine,	//PD line pointer of this DE's Parameter Data Section.
	int ParameterLineCount,	//Line count of this PD lines.
	int DEpointer,			//DEpointer of this DE.
	string* lines	//Output line string variable.Each has the size of 80.
){
	ostringstream sstr1;	
	sstr1<<setw(8)<<m_EntityTypeNumber;//Convert integer to string.
	m_ParameterDataLine=ParameterDataLine;
	sstr1<<setw(8)<<m_ParameterDataLine;
	sstr1<<setw(8)<<MGIges::DEpointer_to_lnumber(m_Structure);
	sstr1<<setw(8)<<m_LineFontPattern;
	sstr1<<setw(8)<<m_Level;
	sstr1<<setw(8)<<MGIges::DEpointer_to_lnumber(m_View);
	sstr1<<setw(8)<<MGIges::DEpointer_to_lnumber(m_TransformationMatrix);
    sstr1<<setw(8)<<MGIges::DEpointer_to_lnumber(m_LabelDisplayAssociativity);
	sstr1<<setw(2)<<setfill('0')<<right<<m_StatusNumber.blankStatus();
	sstr1<<setw(2)<<m_StatusNumber.subordinateEntitySwitch();
	sstr1<<setw(2)<<m_StatusNumber.entityUseFlag();
	sstr1<<setw(2)<<m_StatusNumber.hierarchy();
	sstr1<<'D';
  	int lineN1=MGIges::DEpointer_to_lnumber(DEpointer);	//SequenceNumber output.
 	sstr1<<setw(7)<<setfill('0')<<right<<lineN1;
	string& Oneline=sstr1.str();
	lines[0] = Oneline;

	ostringstream sstrng2;	
    sstrng2<<setw(8)<<m_EntityTypeNumber;
	sstrng2<<setw(8)<<m_LineWeightNumber;
	int ColorNumber=m_ColorNumber;
	if(ColorNumber<0)
		ColorNumber=-DEpointer_to_lnumber(-ColorNumber);
	sstrng2<<setw(8)<<ColorNumber;
	m_ParameterLineCount=ParameterLineCount;
	sstrng2<<setw(8)<<m_ParameterLineCount;
	sstrng2<<setw(8)<<m_FormNumber;
	sstrng2<<setw(8)<<"        ";//m_Reserved1
	sstrng2<<setw(8)<<"        ";//m_Reserved2;
	sstrng2<<setw(8)<<left<<m_EntityLabel;
	sstrng2<<setw(8)<<m_EntitySubscriptNumber;
	sstrng2<<'D';
	int lineN2=lineN1 + 1;
	sstrng2<<setw(7)<<setfill('0')<<right<<lineN2;
	string& Twoline=sstrng2.str();
	lines[1] = Twoline;
}
