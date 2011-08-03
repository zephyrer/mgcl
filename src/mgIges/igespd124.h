/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD124_H__)
#define __MGIGESPD124_H__

#include "mg/Transf.h"
#include "mgiges/IgesPD.h"

///MGIgesPD124 is the class for Iges parameter data type 124(Transformation matrix).
class MGIgesPD124: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD124.
	MGIgesPD124(MGIgesDirectoryEntry* DEpointer=0);
	explicit MGIgesPD124(const MGTransf& tr);

	///Destructor;
	~MGIgesPD124(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

	///convert this transformation to MGTransf.
	void convert_to_MGTransf(MGTransf& tr)const;

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

private:
	///Member data.

	double m_matrix[12];///<(R11,R12,R13,T1,R21,R22,R23,T2,R31,R32,R33,T3)
						///<of the transformation matrix.

	///  | R11 R12 R13 | | X |   | T1 |   | Xout | 
	///  | R21 R22 R23 | | Y | + | T2 | = | Yout |
	///  | R31 R32 R33 | | Z |   | T3 |   | Zout |

};

#endif // __MGIGESPD124_H__
