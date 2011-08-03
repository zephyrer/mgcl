/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD314_H__)
#define __MGIGESPD314_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD314 is the class for Iges parameter data type 314(Color definition entity).
class MGIgesPD314: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD314.
	MGIgesPD314(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD314.
	MGIgesPD314(const MGColor& color);

	///Destructor;
	~MGIgesPD314(){;};

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
		const MGIgesGSec& gsec,	///Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///output plines.
	)const;

public:
//Member data.

	float m_rgb[3];///<RGB percetage data. 0<= m_rgb[.] <=100.
	std::string m_color_name;///<color name.
};

#endif // __MGIGESPD314_H__