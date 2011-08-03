/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD402_H__)
#define __MGIGESPD402_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD402 is the class for Iges parameter data type 402(Group associativity).
class MGIgesPD402: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD402.
	MGIgesPD402(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD402(){;};

	///append a member.
	void append_DE(int de){m_members.push_back(de);};

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

	std::vector<int> m_members;///<vector of directory entries of the group member.
};

#endif // __MGIGESPD402_H__