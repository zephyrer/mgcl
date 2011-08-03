/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __IGES_H__)
#define __IGES_H__

#include <string>
#include <sstream>
#include <ostream>
#include <istream>
#include <iomanip>
#include "mg/Pvector.h"

class MGIgesGSec;
class MGIgesDirectoryEntry;
class MGIgesPD;
class MGIgesParamLine;

///Define namespace for Iges functions.
///all of the number of Iges Parameter Data Section is defined here as a name.
namespace MGIges{

enum EntityTypeNumber{
	Null=0,
	CIRCULAR_ARC=100,
	COMPOSITE_CURVE=102,
	CONIC_ARC=104,
	//CopiousData=106,
	PLANE=108,
	LINE=110,
	PARAMETRIC_SPLINE_CURVE=112,
	//ParametricSplineSurface=114,
	POINT=116,
	RULED_SURFACE=118,
	SURFACE_OF_REVOLUTION=120,
	TABULATED_CYLINDER=122,
	DIRECTION=123,
	TRANSFORMATION_MATRIX=124,
	//Flash=125,
	RATIONAL_BSPLINE_CURVE=126,
	RATIONAL_BSPLINE_SURFACE=128,
	//OffsetCurve=130,
	//ConnectPoint=132,
	//Node=134,
	//FiniteElement=136,
	//NodalDisplacementAndRotation=138,
	//OffsetSurface=140,
	BOUNDARY=141,
	CURVE_ON_PARAMETRIC_SURFACE=142,
	BOUNDED_SURFACE=143,
	TRIMMED_SURFACE=144,
	//NodalResults=146,
	//ElementResults=148,
	//Block=150,
	//RightAngularWedge=152,
	//RightCircularCylinder=154,
	//RightCircularConeFrustum=156,
	SPHERE=158,
	//Torus=160,
	//SolidOfRevolution=162,
	//Ellipsoid=168,
	//BooleanTree=180,
	//SelectedComponent=182,
	//SolidAssembly=184,
	MANIFOLD_SOLID_BREP_OBJECT=186,
	PLANE_SURFACE=190,
	RIGHT_CIRCULAR_CYLINDRICAL_SURFACE=192,
	///RightCircularConicalSurface=194,
	SPHERICAL_SURFACE=196,
	//ToridalSurface=198,
	//AngularDimension=202,
	//CurveDimension=202,
	//CurveDimension=204,
	//DiameterDimension=206,
	//FlagNote=208,
	//GeneralLabel=210,
	//GeneralNote=212,
	//NewGeneralNote=213,
	//Leader=214,
	//LinearDimension=216,
	//OrdinateDimension=218,
	//PointDimension=220,
	//RadiusDimension=222,
	//GeneralSymbol=228,
	//SectionedArea=230,
	//AssociativityDifinition=302,
	//LineFontDefinition=304,
	//MACRODefinition=306,
	//SubfigureDefinition=308,
	//TextFontDefinition=310,
	//TextDisplayTemplate=312,
	COLOR_DEFINITION=314,
	//UnitsData=316,
	//NetworkSubfigureDefinition=320,
	//AttributeTableDefinition=322,
	ASSOCIATIVITY_INSTANCE=402,
	//Drawing=404,
	//Property=406,
	//SingularSubfigureInstance=408,
	//View=410,
	//RectangularArraySubfigureInstance=412,
	//CircularArraySubfigureInstance=414,
	//ExternalReference=416,
	//NodalLoad=418,
	//NetworkSubfigureInstance=420,
	//AttributeTableInstance=422,
	//SolidInstance=430,
	VERTEX=502,
	EDGE=504,
	LOOP=508,
	FACE=510,
	SHELL=514
};

///Convert an MGIgesFstream's DE pointer to the line number to store in IGES file.
///line_number=2*DEpointer-1;
int DEpointer_to_lnumber(int DEpointer);

///Convert a line number stored in IGES file to the MGIgesFstream's DE pointer.
///DEpointer=(line_number+1)/2
int lnumber_to_DEpointer(int line_number);

///Read in DE pointer into DEpointer.
///Line number in the istrm is converted to DE pointer.
///Function's return value is
///  true: when value specified.
///  false:when value not specified, DEpointer be 0.
bool get_DEpointer(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains integer data.
		///<The stream pointer will be advanced to the start position of the next item.
	int& DEpointer	///<output integer data that is converted from the istrm data.
);

///Read in Hollerith_string into strngData.
///Function's return value is
///  true: when value specified, strngData.size() be >0.
///  false:when value not specified, strngData.size() be 0.
bool get_Hollerith_string(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains Hollerith data.
		///<The stream pointer will be advanced to the start position of the next item.
	std::string& strngData///<output string data that is converted from the istrm's Hollerith data.
);

///Read in integer_string into intData.
///Function's return value is
///  true: when value specified.
///  false:when value not specified, intData be 0.
bool get_integer(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains integer data.
		///<The stream pointer will be advanced to the start position of the next item.
	int& intData	///<output integer data that is converted from the istrm data.
);
bool get_integer(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains integer data.
		///<The stream pointer will be advanced to the start position of the next item.
	short& shortData	///<output integer data that is converted from the istrm data.
);

///Read in real_string into realData
///Function's return value is
///  true: when value specified.
///  false:when value not specified, realData be 0.
bool get_real(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains real data.
		///<The stream pointer will be advanced to the start position of the next item.
	double& realData	///<converted real data from istrm will be output.
);
bool get_real(
	char pDelimeter,	///<parameter delimeter
	std::istringstream& istrm,	///<Input string stream that contains real data.
		///<The stream pointer will be advanced to the start position of the next item.
	float& floatData	///<converted real data from istrm will be output.
);

///convert the line id into int(sequence), inputting one line.
void get_ID_sequence(
	const std::string& line,///<Input whole line data(1-80)
	char& sectionID_letter,	///<section identification letter of the line will be output.
	int& sequence			///<ascending sequence number of the line.
);

///Put integer data into plines, converting into string.
///Except for string data, one integer or double data is output
///into one line, not striding over more than one lines.
void put_integer(
	int idata,		///<integer to output.
	const MGIgesGSec& gsec,
	MGPvector<std::string>& plines, ///<output plines.
		///<lines will be added to input plines.
		///<When plines.size()==0, a new-ed std::string whose size()<=63
		///<will be appended. 
		///<When plines.size()>0, let endpline=plines.back(). When (*endpline)
		///<can hold the input data, the data will be appended onto (*endpline).
		///<If one std::string was not enough to store data(greater than 63
		///<characters), new lines of std::string will be created to append after
		///<plines.
	size_t line_len=64///<line length to output, =64(for Parameter data section) or 72.
);

///Put real data into plines, converting into string.
///Except for string data, one integer or double data is output
///into one line, not striding over more than one lines.
void put_real(
	double rdata,	///<double data to output.
	const MGIgesGSec& gsec,
	MGPvector<std::string>& plines, ///<output plines.
		///<lines will be added to input plines.
		///<When plines.size()==0, a new-ed std::string whose size()<=63
		///<will be appended. 
		///<When plines.size()>0, let endpline=plines.back(). When (*endpline)
		///<can hold the input data, the data will be appended onto (*endpline).
		///<If one std::string was not enough to store data(greater than 63
		///<characters), new lines of std::string will be created to append after
		///<plines.
	size_t line_len=64///<line length to output, =64(for Parameter data section) or 72.
);

///Put string data into plines, converting into Hollerith string.
///Only when string data is output(to Holleris string), the data
///may stride over more than one lines.
void put_Hollerith_string(
	const std::string& strngData,///<string data to output,
								///<will be converted to Hollerith data.
	const MGIgesGSec& gsec,
	MGPvector<std::string>& plines, ///<output plines.
		///<lines will be added to input plines.
		///<When plines.size()==0, a new-ed std::string whose size()<=63
		///<will be appended. 
		///<When plines.size()>0, let endpline=plines.back(). When (*endpline)
		///<can hold the input data, the data will be appended onto (*endpline).
		///<If one std::string was not enough to store data(greater than 63
		///<characters), new lines of std::string will be created to append after
		///<plines.
	size_t line_len=64///<line length to output, =64(for Parameter data section) or 72.
);

///Put DE pointer data into plines, converting into line number from DE pointer.
///Except for string data, one integer or double data is output
///into one line, not striding over more than one lines.
///line length to output is always 64(for Parameter data section).
void put_DEpointer(
	int DEpointer,  ///<DE pointer to output.
	const MGIgesGSec& gsec,
	MGPvector<std::string>& plines ///<output plines.
		///<lines will be added to input plines.
		///<When plines.size()==0, a new-ed std::string whose size()<=63
		///<will be appended. 
		///<When plines.size()>0, let endpline=plines.back(). When (*endpline)
		///<can hold the input data, the data will be appended onto (*endpline).
		///<If one std::string was not enough to store data(greater than 63
		///<characters), new lines of std::string will be created to append after
		///<plines.
);

///append record delimeter to plines.
void append_record_delimeter(
		char record_del,
		MGPvector<std::string>& plines
);

};

#endif // __IGES_H__
