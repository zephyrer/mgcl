/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD192_H__)
#define __MGIGESPD192_H__

#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD.h"

///MGIgesPD192 is the class for Iges parameter data type 192
///(Right circular cylindrical surface).
class MGIgesPD192: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD192.
	MGIgesPD192(MGIgesDirectoryEntry* DEpointer=0);

	///Construct PD192, supplying each DE pointer data.
	MGIgesPD192(
		int locationDE, int normalDE, double radius, int refdirDE=0
	);

	///Destructor;
	~MGIgesPD192(){;};

	///Get the plane origin(LOCATION) into origin.
	void getOrigin(const MGIgesIfstream& ifs, MGPosition& origin)const;

	///Get the plane normal into nromal.
	void getNormal(const MGIgesIfstream& ifs, MGUnit_vector& normal)const;

	///Get the plane reference direction(REFDIR) into refdir.
	void getRefdir(const MGIgesIfstream& ifs, MGVector& refdir)const;

	double getRadius()const{return m_radius;};

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

//Member data. These are set as public.

	int m_locationDE;///<a location DE on the cylinder.
	int m_normalDE;///<normal DE of the cylinder(this is a unit vector).
	double m_radius;///<Radius of the cylinder.
	int m_refdirDE;///<reference DE direction of the cylinder if m_refdirDE>0,
				///<=0 if no reference direction.
};

#endif // __MGIGESPD192_H__
