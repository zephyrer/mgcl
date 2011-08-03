/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD104_H__)
#define __MGIGESPD104_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD104 is the class for Iges parameter data type 104(conic arc).
class MGIgesPD104: public MGIgesPD{
public:
	// Constructors.

	/// Constructs an object of class MGIgesPD104.
	MGIgesPD104(MGIgesDirectoryEntry* DEpointer=0);

	///Construct PD100, supplying 2D coordinate data in each array.
	MGIgesPD104(
		const double coef[6],///<Coefficients of the conic equation.
		double Zt,///<Z coordinate of (x,y) plane of the above equation.
		const double start[2], const double terminate[2]
			///<the start and terminate point coordinates of the conic arc in (x,y) plane.
	);

	///Destructor;
	~MGIgesPD104(){;};

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

	double m_coef[6];///<Coefficients of the conic equation.
				///<Let m_coef[]={A,B,C,D,E,F}, then A*x**2+B*x*y+C*y**2+D*x+E*y+F=0.
	double m_Zt;///<Z coordinate of (x,y) plane of the above equation.
	double m_X1, m_Y1;///<the start point coordinate of the conic arc in (x,y) plane.
	double m_X2, m_Y2;///<the terminate point coordinate of the conic arc in (x,y) plane.
};

#endif // __MGIGESPD104_H__