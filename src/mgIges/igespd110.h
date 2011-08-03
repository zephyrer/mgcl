/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD110_H__)
#define __MGIGESPD110_H__

#include "mg/Position.h"
#include "mgiges/IgesPD.h"

class MGIgesDirectoryEntry;

///MGIgesPD110 is the class for Iges parameter data type 110(LINE).
class MGIgesPD110: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD110.
	MGIgesPD110(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD110.
	MGIgesPD110(
		const MGPosition& start,
		const MGPosition& terminate
	);

	///Destructor;
	~MGIgesPD110(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

	///Write out this PD as MGIgesParamLine's(into plines).
	///Except for string data, one integer or double data is output
	///into one MGIgesParamLine, not striding over more than one line.
	///Only when string data is output(to Holleris string), the data
	///may stride over more than one lines.
	///plines[i] for 0<=i<plines.size() are valid.
	void write_out_into_string(
		const MGIgesGSec& gsec,	///<Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///<output plines.
	)const;

public:
//Member data. These are set as public.

	double m_start[3];		///<(x, y, z) of the start point.
	double m_terminate[3];	///<(x, y, z) of the terminate point.
};

#endif // __MGIGESPD110_H__