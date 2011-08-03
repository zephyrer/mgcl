/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD143_H__)
#define __MGIGESPD143_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD143 is the class for Iges parameter data type 143(Bounded Surface).
class MGIgesPD143: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD143.
	MGIgesPD143(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD143(){;};

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
//Member data.

	int m_type;///<The type of bounded surface representation:
			///<=0: the boundary entities shall reference only model space.
			///<=1: The boundary entities shall reference both model space curve
			///<    and the associated parameter space curve collections.
	int m_surface_DE;///<Directory entry of the untrimmed(base) surface.
	std::vector<int> m_boundaries;///<vector of directory entries of the boundary entity(MGIgesPD141).
};

#endif // __MGIGESPD143_H__