/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD112_H__)
#define __MGIGESPD112_H__

#include <vector>
#include "mg/NDDArray.h"
#include "mgiges/IgesPD.h"

///@cond
//A Private class for MGIgesPD112, express one coefficient of IGES NURBS.
class MGIgesSplineCoef{
public:
	double m_abcd[4][3];
	///<coefficients of one segment=m_abcd[.][i]=(A,B,C,D) for i-th space dimension, that is,
	///<m_abcd[.][0]=(A,B,C,D) for X, m_abcd[.][1] for Y, and m_abcd[.][2] for Z.
};
///@endcond

///MGIgesPD112 is the class for Iges parameter data type 112(Parametric spline curve).
class MGIgesPD112: public MGIgesPD{
public:
	// Constructors.

	/// Constructs an object of class MGIgesPD112.
	MGIgesPD112(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD112(){;};

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

	short m_spline_type;	///<=1:linear, =2:quardric, =3:cubic,
						///<=4:Wilson-Fowler, =5:Modified Wilson-Fowler, =6:B-spline.
	short m_continuity;	///<Degree of continuity with respect to arc length.
	short m_dimension;	///<=2:planar, =3:nonplanar.
	MGNDDArray m_tau;
				///<Break point sequence of the spline. m_tau.length()=number_of_segments+1;
	std::vector<MGIgesSplineCoef> m_coefs;
				///<Coefficients sequence of the spline. m_coefs.size()=number_of_segments+1;
			
};

#endif // __MGIGESPD112_H__