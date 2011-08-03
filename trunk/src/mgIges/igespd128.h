/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD128_H__)
#define __MGIGESPD128_H__

#include <vector>
#include "mg/Position.h"
#include "mg/KnotVector.h"
#include "mg/SPointSeq.h"
#include "mgiges/IgesPD.h"

///MGIgesPD128 is the class for Iges parameter data type 128(NURBS Surface).
class MGIgesPD128: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD128.
	MGIgesPD128(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD128.
	MGIgesPD128(const MGSBRep& sb);

	/// Constructs an object of class MGIgesPD128.
	MGIgesPD128(const MGRSBRep& sb);

	///Destructor;
	~MGIgesPD128(){;};

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

///Member data.

	///Here we denote the surface S(u,v) using the parameter (u,v).
	///nBrepU=m_upper_indexU+1, orderU=m_degreeU+1
	///nBrepV=m_upper_indexV+1, orderV=m_degreeV+1
	///(nBrep is B-Representation dimension of B-Spline).

	int m_upper_indexU;///<Upper index of sum along the parameter u,
						///<that is, m_upper_indexU=nBrepU-1.
	int m_upper_indexV;///<Upper index of sum along the parameter v,
						///<that is, m_upper_indexV=nBrepV-1.
	int m_degreeU;	///<Degree of the NURBS along the parameter u,
					///<that is, m_degreeU=orderU-1.
	int m_degreeV;	///<Degree of the NURBS along the parameter v,
					///<that is, m_degreeV=orderV-1.
	short m_closedU;///<=0:open , =1:closed curve, along the parameter u;	
	short m_closedV;///<=0:open , =1:closed curve, along the parameter v;	
	short m_non_rational;///<=0:rational, =1:non rational;	
	short m_periodicU;	///<=0:nonperiodic, =1:periodic, along the parameter u;
	short m_periodicV;	///<=0:nonperiodic, =1:periodic, along the parameter v;
	MGKnotVector m_knotsU;///<Knot vector of length (nBrepU+orderU).
	MGKnotVector m_knotsV;///<Knot vector of length (nBrepV+orderV).
	MGSPointSeq m_weights;///<Weight vector of length nBrepU*nBrepV.
		///<m_weights[i+nBrepU*j] is the weight for the control points m_control_points[i+nBrepU*j],
		///<for 0<=i<=nBrepU-1, 0<=j<=nBrepV-1.
	MGSPointSeq m_control_points;///<Control points of length nBrepU*nBrepV.
		///<m_control_points[i+nBrepU*j] for 0<=i<=nBrepU-1, 0<=j<=nBrepV-1.
	double m_start_paramU, m_end_paramU;///<Starting and ending parameters along the parameter U,
		///<that is, m_start_paramU=m_knotsU[m_degreeU], or m_start_paramU=m_knots[orderU-1].
		///<m_end_paramU=m_knots[nBrepU],
	double m_start_paramV, m_end_paramV;///<Starting and ending parameters along the parameter v,
		///<that is, m_start_paramV=m_knotsV[m_degreeV], or m_start_paramV=m_knots[orderV-1].
		///<m_end_paramV=m_knots[nBrepV],
};

#endif // __MGIGESPD128_H__