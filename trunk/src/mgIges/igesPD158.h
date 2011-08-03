/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD158_H__)
#define __MGIGESPD158_H__

#include <vector>
#include "mg/Sphere.h"
#include "mgiges/IgesPD.h"

///MGIgesPD158 is the class for Iges parameter data type 158(unparameterised sphere).
class MGIgesPD158: public MGIgesPD{
public:
	// Constructors.

	/// Constructs an object of class MGIgesPD158.
	MGIgesPD158(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD158 from an MGSphere
	MGIgesPD158(const MGSphere& sphere);

	///Destructor;
	~MGIgesPD158(){;};

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

	double m_radius;	///<radius of the sphere.
	double m_center[3];///<center coordinates of the sphere.
};

#endif // __MGIGESPD158_H__