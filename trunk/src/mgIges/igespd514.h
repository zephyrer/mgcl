/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD514_H__)
#define __MGIGESPD514_H__

#include <vector>
#include "mg/Position.h"
#include "mgiges/IgesPD.h"

///MGIgesPD514 is the class for Iges parameter data type 514(Shell).
class MGIgesPD514: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD514.
	MGIgesPD514(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD514(){;};

	///append a face.
	void push_back(
		int face_DE,	///DE pointer of the face.
		bool same_direction=true
	);

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
	bool m_is_closed;///<Indicated if this is an closed shell(true) or open(false).
	std::vector<int> m_faces;  ///<vector of the face(IgesPD510) DE pointers.
	std::vector<bool> m_orientations;///<vector of the bools that indicate
		///<whether m_face[i]'s direction agrees with the direction of
		///<the underlying surface(=true) or not.
};

#endif // __MGIGESPD514_H__
