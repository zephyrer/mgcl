/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD123_H__)
#define __MGIGESPD123_H__

#include "mgiges/IgesPD.h"

class MGVector;

///MGIgesPD123 is the class for Iges parameter data type 123(DIRECTION).
class MGIgesPD123: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD123.
	MGIgesPD123(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD123.
	MGIgesPD123(const MGVector& vec);

	///Destructor;
	~MGIgesPD123(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

	///Convert the direction data to MGVector direction.
	void convert_to_vector(MGVector& direction)const;

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

private:
	//Member data. These are set as public.

	double m_xyz[3];///<(x,y,z) coordinates of the direction vector.
};

#endif // __MGIGESPD123_H__