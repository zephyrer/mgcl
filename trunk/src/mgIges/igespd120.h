/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD120_H__)
#define __MGIGESPD120_H__

#include "mgiges/IgesPD.h"

///MGIgesPD120 is the class for Iges parameter data type 120(Surface of Revolution).
class MGIgesPD120: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD120.
	MGIgesPD120(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD120(){;};

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

	int m_axis_of_revolution_DE;///<Directory entry of the line entity of axis of revolution.
	int m_generatrix_DE;		///<Directory entry of the generatrix.
	double m_start_angle;///<start angle of revolution in radian.
	double m_terminate_angle;///<terminate angle of revolution in radian.
};

#endif // __MGIGESPD120_H__