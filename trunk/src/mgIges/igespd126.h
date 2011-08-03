/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD126_H__)
#define __MGIGESPD126_H__

#include "mg/Position.h"
#include "mg/LBRep.h"
#include "mgiges/IgesPD.h"

///MGIgesPD126 is the class for Iges parameter data type 126(NURBS).
class MGIgesPD126: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD126.
	MGIgesPD126(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD126.
	MGIgesPD126(const MGLBRep& lb);

	/// Constructs an object of class MGIgesPD126.
	MGIgesPD126(const MGRLBRep& lb);

	///Destructor;
	~MGIgesPD126(){;};

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

//Member data.

	///Here we denote nBrep=m_upper_index+1, order=m_degree+1
	///(nBrep is B-Representation dimension of B-Spline).

	int m_upper_index;///<Upper index of sum, that is, m_upper_index=nBrep-1.
	int m_degree;	///<Degree of the NURBS, that is, m_degree=order-1.
	short m_planar;	///<=0:nonplanar, =1:planar;	
	short m_closed;	///<=0:open curve, =1:closed curve;	
	short m_non_rational;///<=0:rational, =1:non rational;	
	short m_periodic;	///<=0:nonperiodic, =1:periodic;
	MGKnotVector m_knots;///<Knot vector of length (nBrep+order).
	std::vector<double> m_weights;///<Weight vector of length nBrep.
	MGBPointSeq m_control_points;///<Control points of length nBrep.
	double m_start_param, m_end_param;///<Starting and ending parameters, that is,
			///<m_start_param=m_knots[m_degree], or m_start_param=m_knots[order-1].
			///<m_end_param=m_knots[nBrep],
	double m_normal[3];///<Normal vector of the plane if the NURBS is planar.
};

#endif // __MGIGESPD126_H__
