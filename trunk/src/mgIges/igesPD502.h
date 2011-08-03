/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD502_H__)
#define __MGIGESPD502_H__

/// @file
///	@brief  Declaration for class MGIgesPD502.
///	@author DG Technologies(http:///www.dgtech.co.jp/)

#include <vector>
#include "mg/Position.h"
#include "mgiges/IgesPD.h"

///MGIgesPD502 is the class for
///the Iges parameter data type 502(VERTEX List Entity) form 1.
class MGIgesPD502: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD502.
	MGIgesPD502(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD502(){;};

	MGPosition& operator[](int i){return m_vertices[i];};
	const MGPosition& operator[](int i)const{return m_vertices[i];};

	///append one vertex data.
	void push_back(const MGPosition& vertex);

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

	///Vertices of 3D coordinates.
	std::vector<MGPosition> m_vertices;///<m_vertices[0] is dummy. This is because
			///<list index of MGIges504Edge's m_Svertex/m_Tvertex starts from 1.
};

#endif // __MGIGESPD502_H__
