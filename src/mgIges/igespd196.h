/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD196_H__)
#define __MGIGESPD196_H__

#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Sphere.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD.h"

///MGIgesPD196 is the class for Iges parameter data type 196(sphere surface)
class MGIgesPD196: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD196.
	MGIgesPD196(MGIgesDirectoryEntry* DEpointer=0);

	///Construct PD196, supplying each DE pointer data.
	MGIgesPD196(
		int locationDE, double radius, int axisDE=0, int refdirDE=0
	);

	///Destructor;
	~MGIgesPD196(){;};

	///Get the sphere center(LOCATION) into origin.
	void getCenter(const MGIgesIfstream& ifs, MGPosition& center)const;

	///Get the plane normal into nromal.
	void getAxis(const MGIgesIfstream& ifs, MGUnit_vector& axis)const;

	///Get the sphere reference direction(REFDIR) into refdir.
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
		const MGIgesGSec& gsec,	///<Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///<output plines.
	)const;

	///Member data. These are set as public.

	int m_locationDE;///<center DE of the sphere
	double m_radius;///<Radius of the sphere.
	int m_axisDE;///<normal DE of the sphere(this is a unit vector) if m_axisDE>0,
				///<=0 if no axis direction.
	int m_refdirDE;///<reference DE direction of the sphere if m_refdirDE>0,
				///<=0 if no reference direction.
};

#endif // __MGIGESPD196_H__
