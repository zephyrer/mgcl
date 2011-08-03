/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD122_H__)
#define __MGIGESPD122_H__

#include "mgiges/IgesPD.h"

///MGIgesPD122 is the class for Iges parameter data type 122(Tabulated Cylinder).
class MGIgesPD122: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD122.
	MGIgesPD122(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD122.
	MGIgesPD122(int diretrix_DE, double terminate_point[3]);

	///Destructor;
	~MGIgesPD122(){;};

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

	int m_directrix_DE;///<Directory entry of the directrix curve.
	double m_terminate_point[3];///<(x,y,z) coordinates of the terminate point of
						///<the start point of the generatrix.
};

#endif // __MGIGESPD122_H__