/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD118_H__)
#define __MGIGESPD118_H__

#include "mgiges/IgesPD.h"

///MGIgesPD118 is the class for Iges parameter data type 118(Ruled Surface).
class MGIgesPD118: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD118.
	MGIgesPD118(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD118(){;};

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

	int m_1st_Curve_DE;	///<Directory entry of the 1st curve of the ruled surface.
	int m_2nd_Curve_DE;	///<Directory entry of the 2nd curve of the ruled surface.
	short m_direction_flag;///<direction flag,
						///<=0:Join 1st to 1st and last to last,
						///<=1:Join 1st to last, and last to 1st.
	short m_developable_flag;///<=0:Possibly not developable, =1:developable.
};

#endif // __MGIGESPD118_H__
