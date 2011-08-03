/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESGSEC_H__)
#define __MGIGESGSEC_H__

#include <string>
class MGIgesOfstream;

/** @addtogroup FileInputOutput
 *  @{
 */

///MGIgesGSec describes a Global Section of a IGES file.
class MGIgesGSec{
/// Constructors.

public:
	/// Constructs an object of class MGIgesGSec.
	///Default constructor, includes all the defalut value of MGCL.
	MGIgesGSec(
		const char* filename=0	///<Input IGES file name
	);

	/// Read MGIgesGSec data into this object from a global section
	///string that includes all the string of a IGES Global section.
	///gsec_string includes all the global section string that does not
	///have the identification codes and the sequence number.
	///each item in gsec_string is separated by the parameter delimeter charactors,
	///and the end character of gsec_string is the record delimeter  charactor.
	void read_in(const std::string& gsec_string);

	///Write out this Global section to MGIgesOfstream.
	///Return is the number of lines output.
	int write_out(MGIgesOfstream& ofs);

	char paramDelimeter()const{return m_delimeter_param;};
	char recordDelimeter()const{return m_delimeter_record;};

////Member data. These are set as public.
	char m_delimeter_param;						// 1
	char m_delimeter_record;					// 2
	std::string m_productID_sender;				// 3
	std::string m_file_name;					// 4
	std::string m_native_systemID;				// 5
	std::string m_preprocessor_version;			// 6
	int m_number_of_bits_of_integer;			// 7
	int m_magnitude_single_precision;			// 8
	int m_significance_single_precision;		// 9
	int m_magnitude_double_precision;			//10
	int m_significance_double_precision;		//11
	std::string m_productID_receiver;			//12
	double m_model_space_scale;					//13
	int m_unit_flag;							//14
	std::string m_unit_name;					//15
	int m_max_number_of_line_weight_gradations;	//16
	double m_width_of_max_line_weight;			//17
	std::string m_DateTime_File_generation;		//18
	double m_min_resolution;					//19
	double m_max_coordinate_value;				//20
	std::string m_author_name;					//21
	std::string m_author_organazation;			//22
	int m_version_flag;							//23
	int m_drafting_standard_flag;				//24
	std::string m_DateTime_Model_generation;	//25
	std::string m_application_protocolID;		//26
};

/** @} */ // end of FileInputOutput group
#endif // __MGIGESGSEC_H__
