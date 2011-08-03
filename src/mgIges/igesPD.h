/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD_H__)
#define __MGIGESPD_H__

#include "mgiges/Iges.h"

// forward declerations
class MGIgesDirectoryEntry;
class MGIgesParamLine;

///MGIgesPD is the parent class of all the Parameter data section type.
///Each type of parameter data section will be inheritted from this class.
class MGIgesPD{
	friend MGIgesDirectoryEntry;

/// Constructors.
public:
	/// Constructs an object of class MGIgesPD.
	MGIgesPD();///Default constructor.

	/// Constructs an object of class MGIgesPD.
	MGIgesPD(int type_number, MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	virtual ~MGIgesPD();

	///Read in parameter data from string stream data.
	virtual void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	)=0;

	void setDE(MGIgesDirectoryEntry* DE){m_DEpointer=DE;};
	int type_number()const{return m_type_number;};
	const MGIgesDirectoryEntry* DEpointer()const{return m_DEpointer;};
	MGIgesDirectoryEntry* DEpointer(){return m_DEpointer;};

	///Write out this PD as MGIgesParamLine's(into plines).
	///Except for string data, one integer or double data is output
	///into one MGIgesParamLine, not striding over more than one line.
	///Only when string data is output(to Holleris string), the data
	///may stride over more than one lines.
	///plines[i] for 0<=i<plines.size() are valid.
	virtual void write_out_into_string(
		const MGIgesGSec& gsec,	///<Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///<output plines.
	)const=0;

private:
//Member data. These are set as public.

	int m_type_number;
	MGIgesDirectoryEntry* m_DEpointer;	///<DE pointer of this parameter data.
};

#endif // __MGIGESPD_H__
