/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD100_H__)
#define __MGIGESPD100_H__

#include "mgiges/IgesPD.h"

///MGIgesPD100 is the class for Iges parameter data type 100(circular arc).
class MGIgesPD100: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD100.
	MGIgesPD100(MGIgesDirectoryEntry* DEpointer=0);

	///Construct PD100, supplying 2D coordinate data in each array.
	MGIgesPD100(
		const double center[2], const double start[2], const double terminate[2],
		double Z=0.	///Z coord
	);

	///Destructor;
	~MGIgesPD100(){;};

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

//Member data. These are set as public.
	double m_Zt;
	double m_center[2];///<(x1, y1) of the center.
	double m_start[2];///<(x2, y2) of the start point.
	double m_terminate[2];///<(x2, y2) of the terminate point.
};

#endif // __MGIGESPD100_H__
