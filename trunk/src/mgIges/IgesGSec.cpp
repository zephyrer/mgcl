/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesGSec.
//!	@author DG Technologies(http://www.dgtech.co.jp/)
#include "MGCLStdAfx.h"
#include "mgGL/LineWidth.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//!	@brief MGIgesGSec describes a Global Section of a IGES file.
using namespace MGIges;

//! Constructs an object of class MGIgesGSec.
//Default constructor, includes all the defalut value of MGCL.
MGIgesGSec::MGIgesGSec(
	const char* filename	//Input IGES file name
):m_delimeter_param(','), m_delimeter_record(';'), m_productID_sender("MGCL"),
m_native_systemID("MGCL Version 7.1"),
m_preprocessor_version("IGES 5.3"), m_number_of_bits_of_integer(32),
m_magnitude_single_precision(34), m_significance_single_precision(7),
m_magnitude_double_precision(34), m_significance_double_precision(15),
m_model_space_scale(1.), m_unit_flag(2), m_unit_name("MM"),
m_min_resolution(.001), m_max_coordinate_value(99000.),m_author_name("  "),
m_author_organazation("DG Technologies, Inc."),m_version_flag(701),
m_drafting_standard_flag(0){
	if(filename)
		m_file_name=std::string(filename);
	else
		m_file_name=std::string("MGCLIGES.iges");

	const time_t timer = ::time(0);
	const tm* date = ::localtime(&timer);
	const size_t nbuf = 256;
	char fmtdat[nbuf] = {0};
	const size_t nlen = ::strftime(fmtdat, nbuf - 1, "%Y%m%d.%H%M%S", date);

	m_DateTime_Model_generation.assign(&fmtdat[0], nlen);
	m_DateTime_File_generation.assign(&fmtdat[0], nlen);

	MGLineWidth w;
	float fwidth=w.get_maximum_width();
	m_max_number_of_line_weight_gradations=int(fwidth);
	m_width_of_max_line_weight=m_max_number_of_line_weight_gradations;
}

//! Read MGIgesGSec data into this object from a global section
//string that includes all the string of a IGES Global section.
//gsec_string includes all the global section string that does not
//have the identification codes and the sequence number.
//each item in gsec_string is separated by the parameter delimeter charactors,
//and the end character of gsec_string is the record delimeter  charactor.
void MGIgesGSec::read_in(const std::string& gsec_string){
	std::istringstream gsstream(gsec_string);
	std::string onep;

	//1. parameter delimeter.
	char pdelim;
	gsstream.read(&pdelim,1);
	if(pdelim != ','){
		assert(pdelim=='1');
		char H;
		gsstream.read(&H,1); assert(H=='H');
		gsstream.read(&m_delimeter_param,1);
		gsstream.read(&pdelim,1);assert(pdelim==m_delimeter_param);
	}

	//2. record delimeter.
	get_Hollerith_string(m_delimeter_param,gsstream,onep);
	if(onep.length()){
		assert(onep.length()==1);
		m_delimeter_record=onep[0];
	}else
		m_delimeter_record=';';

	//3. productID_sender
	get_Hollerith_string(m_delimeter_param,gsstream,m_productID_sender);

	//4. file_name
	get_Hollerith_string(m_delimeter_param,gsstream,m_file_name);

	//5. native_systemID
	get_Hollerith_string(m_delimeter_param,gsstream,m_native_systemID);

	//6. m_preprocessor_version
	get_Hollerith_string(m_delimeter_param,gsstream,m_preprocessor_version);

	//7. m_number_of_bits_of_integer
	get_integer(m_delimeter_param,gsstream,m_number_of_bits_of_integer);

	//8. m_magnitude_single_precision
	get_integer(m_delimeter_param,gsstream,m_magnitude_single_precision);

	//9. m_significance_single_precision
	get_integer(m_delimeter_param,gsstream,m_significance_single_precision);

	//10. m_magnitude_double_precision
	get_integer(m_delimeter_param,gsstream,m_magnitude_double_precision);

	//11. m_significance_double_precision
	get_integer(m_delimeter_param,gsstream,m_significance_double_precision);

	//12. m_productID_receiver
	if(!get_Hollerith_string(m_delimeter_param,gsstream,m_productID_receiver))
		m_productID_receiver=m_productID_sender;

	//13. m_model_space_scale
	if(!get_real(m_delimeter_param,gsstream,m_model_space_scale))
		m_model_space_scale=1.;

	//14. m_unit_flag
	if(!get_integer(m_delimeter_param,gsstream,m_unit_flag))
		m_unit_flag=1;

	//15. m_unit_name
	if(!get_Hollerith_string(m_delimeter_param,gsstream,m_unit_name))
		m_unit_name="INCH";

	//16. m_max_number_of_line_weight_gradations
	if(!get_integer(m_delimeter_param,gsstream,m_max_number_of_line_weight_gradations))
		m_max_number_of_line_weight_gradations=1;

	//17. m_width_of_max_line_weight
	get_real(m_delimeter_param,gsstream,m_width_of_max_line_weight);

	//18. m_DateTime_File_generation
	get_Hollerith_string(m_delimeter_param,gsstream,m_DateTime_File_generation);

	//19. m_min_resolution
	get_real(m_delimeter_param,gsstream,m_min_resolution);

	//20. m_max_coordinate_value
	get_real(m_delimeter_param,gsstream,m_max_coordinate_value);

	//21. m_author_name
	get_Hollerith_string(m_delimeter_param,gsstream,m_author_name);

	//22. m_author_organazation
	get_Hollerith_string(m_delimeter_param,gsstream,m_author_organazation);

	//23. m_version_flag
	get_integer(m_delimeter_param,gsstream,m_version_flag);
	if(m_version_flag<1)
		m_version_flag=3;
	else if(m_version_flag>11)
		m_version_flag=11;

	//24. m_drafting_standard_flag
	get_integer(m_delimeter_param,gsstream,m_drafting_standard_flag);

	//25. m_DateTime_Model_generation
	get_Hollerith_string(m_delimeter_param,gsstream,m_DateTime_Model_generation);

	//26. m_application_protocolID
	get_Hollerith_string(m_delimeter_param,gsstream,m_application_protocolID);
}
