/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD190_H__)
#define __MGIGESPD190_H__

#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD.h"

///MGIgesPD190 is the class for Iges parameter data type 190(plane surface).
class MGIgesPD190: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD190.
	MGIgesPD190(MGIgesDirectoryEntry* DEpointer=0);

	///Construct PD190, supplying each DE pointer data.
	MGIgesPD190(
		int locationDE, int normalDE, int refdirDE=0
	);

	///Destructor;
	~MGIgesPD190(){;};

	///Get the plane origin(LOCATION) into origin.
	void getOrigin(const MGIgesIfstream& ifs, MGPosition& origin)const;

	///Get the plane normal into nromal.
	void getNormal(const MGIgesIfstream& ifs, MGUnit_vector& normal)const;

	///Get the plane reference direction(REFDIR) into refdir.
	void getRefdir(const MGIgesIfstream& ifs, MGVector& refdir)const;

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
	int m_locationDE;///<a location DE on the plane.
	int m_normalDE;///<normal DE of the plane(this is a unit vector).
	int m_refdirDE;///<reference DE direction of the plane if m_refdirDE>0,
				///<=0 if no reference direction.
};

#endif // __MGIGESPD190_H__
