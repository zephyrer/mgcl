/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD102_H__)
#define __MGIGESPD102_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD102 is the class for Iges parameter data type 102(Composite curve).
class MGIgesPD102: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD102.
	MGIgesPD102(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD102(){;};

	///append a new curve.
	void append_curve(int curve_de){m_curve_DEs.push_back(curve_de);};

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

	std::vector<int> m_curve_DEs;	///<pointer to each member curve's directory entry.
};

#endif // __MGIGESPD102_H__