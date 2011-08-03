/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGEDIRECTORYENTRY_H__)
#define __MGIGEDIRECTORYENTRY_H__

#include <memory>
#include <string>
#include "mgiges/IgesGSec.h"
#include "mgiges/IgesPD.h"
#include "mgiges/IgesDEStatusNumber.h"

MGEXTERN const MGIgesDEStatusNumber mgIgesDEStatusNumber;

///	@brief MGIgesDirectoryEntry describes a directory entry section of an IGES file.
///In IGES file one directory entry is stored in two records. MGIgesDirectoryEntry
///holds the two records data.
class MGIgesDirectoryEntry{
/// Constructors.
public:
	/// Constructs an object of class MGIgesDirectoryEntry.
	MGIgesDirectoryEntry();///Default constructor, includes all the defalut value
				///of MGCL.
	MGIgesDirectoryEntry(
		int EntityTypeNumber,
		const std::string& EntityLabel,
		const MGIgesDEStatusNumber& StatusNumber=mgIgesDEStatusNumber,
		int ColorNumber=0,///<When negated value, is a pointer to the directory entry.
		int LineWeightNumber=0,
		int LineFontPattern=0,///<When negated value, is a pointer to the directory entry.
		int TRANSFORMATION_MATRIX=0,///<Pointer to the directory entity of the matrix.
		int FormNumber=0,
		int EntitySubscriptNumber=0
	);

	///Construct from a string that contais one pair of lines of IGES File.
	MGIgesDirectoryEntry(const std::string& DEstring);

	int EntityTypeNumber()const{return m_EntityTypeNumber;};
	bool is_visible()const{return m_StatusNumber.blankStatus()==0;};
	bool is_independent()const{return m_StatusNumber.subordinateEntitySwitch()==0;};
	const MGIgesDEStatusNumber& status_number()const{return m_StatusNumber;};
	int ParameterDataLine()const{return m_ParameterDataLine;};
	int ParameterLineCount()const{return m_ParameterLineCount;};
	const std::auto_ptr<MGIgesPD>& paramData()const{return m_ParamData;};
	int transformID()const{return m_TransformationMatrix;};
	int LineFontPattern()const{return m_LineFontPattern;};
	int LineWeightNumber()const{return m_LineWeightNumber;};
	float LineWidth(const MGIgesGSec& gsec)const;
	int ColorNumber()const{return m_ColorNumber;};
	int FormNumber()const{return m_FormNumber;};

	void setFormNumber(int FormNumber){m_FormNumber=FormNumber;};
	void setPD(std::auto_ptr<MGIgesPD>& pd);
	void set_SubordinateEntitySwitch(MGIgesDEStatusNumber::SESwitch eswitch){
		m_StatusNumber.set_SubordinateEntitySwitch(eswitch);
	};
	void setTransformID(int trid){m_TransformationMatrix=trid;};
	void setLineWeightNumber(int weight){m_LineWeightNumber=weight;};

	///Convert paramData to the each type parameter data, and set it in m_ParamData.
	void setPD(
		char pDelimeter,		///<parameter delimeter
		const std::string& paramData	///<whole parameter data string,
							///<Each of the parameter is separated by pDelimeter.
	);

	///Put this DE data to 2 lines of IGES file string data.
	void put_to_string(
		int ParameterDataLine,	///<PD line pointer of this DE's Parameter Data Section.
		int ParameterLineCount,	///<Line count of this PD lines.
		int DEpointer,			///<DEpointer of this DE.
		std::string lines[2]	///<Output line string variable,Each has the size of 80.
	);

private:
	///Set pd in m_ParamData.
	void setParamData(std::auto_ptr<MGIgesPD>& pd=std::auto_ptr<MGIgesPD>(0));

private:
	int m_EntityTypeNumber;
	int m_ParameterDataLine;///<Poiter to the 1st line of the parameter data record.
	int m_Structure;	///<When negated value, is a pointer to the directory entry.
	int m_LineFontPattern;///<When negated value, is a pointer to the directory entry.
	int m_Level;///<When negated value, is a pointer to the directory entry.
	int m_View;///<Pointer to the directory entity.
	int m_TransformationMatrix;///<directory entity Pointer of the transform.
	int m_LabelDisplayAssociativity;
	MGIgesDEStatusNumber m_StatusNumber;
	int m_LineWeightNumber;
	int m_ColorNumber;///<When negated value, is a pointer to the directory entry.
	int m_ParameterLineCount;///<Line count of Prameter Data Lines.
	int m_FormNumber;
	//int m_Reserved1;
	//int m_Reserved2;
	std::string m_EntityLabel;
	int m_EntitySubscriptNumber;

	std::auto_ptr<MGIgesPD> m_ParamData;
};

#endif // __MGIGEDIRECTORYENTRY_H__
