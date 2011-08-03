/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD116_H__)
#define __MGIGESPD116_H__

#include "mg/Point.h"
#include "mgiges/IgesPD.h"

///MGIgesPD116 is the class for Iges parameter data type 116(POINT).
class MGIgesPD116: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD116.
	MGIgesPD116(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD116.
	MGIgesPD116(const MGPoint& P,int display_symbolDE=0);
	MGIgesPD116(const MGPosition& P,int display_symbolDE=0);

	/// Constructs an object of class MGIgesPD116.
	MGIgesPD116(const double coordinates[3],int display_symbolDE=0);

	///Destructor;
	~MGIgesPD116(){;};

	///Convert the point data to MGPosition position.
	void convert_to_position(MGPosition& position)const;

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

	double m_coordinates[3];///<(x, y, z) of the point.
	int m_display_symbolDE;	///<Directory entry of the subfigure definition of the display symbol.
};

#endif // __MGIGESPD116_H__