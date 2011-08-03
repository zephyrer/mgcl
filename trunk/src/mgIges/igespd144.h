/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD144_H__)
#define __MGIGESPD144_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD144 is the class for Iges parameter data type 144(Trimmed Surface).
class MGIgesPD144: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD144.
	MGIgesPD144(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD144.
	MGIgesPD144(
		int surfaceDE,		///<Base surface DE.
		int outerboundaryDE	///<if =0, no outer boundary.
	);

	///Destructor;
	~MGIgesPD144(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

	void append_inner_boundary(int inner_boundaryDE){
		m_inner_boundaries.push_back(inner_boundaryDE);
	};

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
//Member data.

	int m_surface_DE;///<Directory entry of the untrimmed(base) surface.
	int m_outer_boundary_type;
					///<=0: outer boundary is the boudary of m_surface_DE.
					///<=1: otherwise.
	int m_outer_boudary_DE;///<outer boundary DE of the parametric space curve.
	std::vector<int> m_inner_boundaries;///<vector of directory entry of the inner boundary entities.
};

#endif // __MGIGESPD144_H__