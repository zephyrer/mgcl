/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD108_H__)
#define __MGIGESPD108_H__

#include <vector>
#include "mg/Plane.h"
#include "mgiges/IgesPD.h"

///MGIgesPD108 is the class for Iges parameter data type 108(Plane).
class MGIgesPD108: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD108.
	MGIgesPD108(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD108 from a MGPlane
	MGIgesPD108(const MGPlane& plane);

	///Destructor;
	~MGIgesPD108(){;};

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
///Member data. These are set as public.

	double m_coef[4];///<Plane coefficients m_coef[.]={A,B,C,D} where A*X+B*Y+C*Z=D;
	int m_boundCurve_DE;///<directory entry of bounding curve, maybe null.
	double m_ref_point[3];///<Reference point on the plane(at which symbol be displayed).
	double m_symbol_size;///<Symbol size to display.
};

#endif // __MGIGESPD108_H__