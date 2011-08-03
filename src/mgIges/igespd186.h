/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD186_H__)
#define __MGIGESPD186_H__

#include <vector>
#include "mgiges/IgesPD.h"

///MGIgesPD186 is the class for Iges parameter data type 186
///(MSBO:Manifold Solid B-Rep Object Entity).
class MGIgesPD186: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD186.
	MGIgesPD186(MGIgesDirectoryEntry* DEpointer=0)
		:MGIgesPD(MGIges::MANIFOLD_SOLID_BREP_OBJECT,DEpointer){;};

	/// Constructs an object of class MGIgesPD186.
	MGIgesPD186(
		int shellDE,	///SHELL DE.
		bool orientation=true///Orientation flag of shell with repsect to
				///its underlying faces, =true:agrees
	):m_shell_DE(shellDE),m_orientation(orientation){;};

	///Destructor;
	~MGIgesPD186(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

	void append_void_shell(
		int void_shell_DE,
		bool orientation=true
	){
		m_void_shells.push_back(void_shell_DE);
		m_orientations.push_back(orientation);
	};

	///Write out this PD as MGIgesParamLine's(into plines).
	///Except for string data, one integer or double data is output
	///into one MGIgesParamLine, not striding over more than one line.
	///Only when string data is output(to Holleris string), the data
	///may stride over more than one lines.
	///plines[i] for 0<=i<plines.size() are valid.
	void write_out_into_string(
		const MGIgesGSec& gsec,	///Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///output plines.
	)const;

public:
//Member data.
	int m_shell_DE;///<Directory entry of the untrimmed(base) surface.
	bool m_orientation;
			///<=true: the shell orientation agrees to its underlying faces.
	std::vector<int> m_void_shells;///<vector of void shells.
	std::vector<bool> m_orientations;
			///<m_orientations[i] is the orientaion of the i-th void shell m_void_shells[i].
};

#endif // __MGIGESPD186_H__