/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD510_H__)
#define __MGIGESPD510_H__

#include <vector>
#include "mg/Position.h"
#include "mgiges/IgesPD.h"
class MGIgesIfstream;
class MGIgesOfstream;
class MGFace;

///MGIgesPD510 is the class for Iges parameter data type 510(FACE).
class MGIgesPD510: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD510.
	MGIgesPD510(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD510(){;};

	///consvert m_surface_DE surface to MGSurface.
	///Returned is a newed MGSurface object.
	MGSurface* convert_to_surface(const MGIgesIfstream& igesIstream)const;

	///append an edge.
	void push_back(int loop){m_loops.push_back(loop);};

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

	int m_surface_DE;///<pointer to the DE of the underlying surface.
	bool m_outer_loop_identified;///<Indicates if outer loop is identified(true) or not.
			///<When true, 
	std::vector<int> m_loops;///<Pointers to the DE that constitue the boundary of this face.
			///<When m_outer_loop_identified=true, m_loops[0] is the outer loop, else
			///<all of the loops are inner loops.
};

#endif // __MGIGESPD510_H__
