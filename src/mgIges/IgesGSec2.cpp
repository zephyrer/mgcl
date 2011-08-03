/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesGSec.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesOfstream.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
using namespace MGIges;

//!	@brief MGIgesGSec describes a Global Section of a IGES file.

//Write out this Global section to MGIgesOfstream.
//Return is the number of lines output.
int MGIgesGSec::write_out(MGIgesOfstream& ofs){
		
	int lineID=1;
	std::ofstream& ofsr=ofs.ofstrm();//Output stream.
	MGPvector<std::string> plines;

	//1.parameter delimeter.
	std::string delp(1,m_delimeter_param);
	put_Hollerith_string(delp,*this,plines,72);
	
	//2.record delimeter.
	std::string delr(1,m_delimeter_record);
	put_Hollerith_string(delr,*this,plines,72);	
	
	//3.productID_ender.
	put_Hollerith_string(m_productID_sender,*this,plines,72);
	
	//4.file_name
	put_Hollerith_string(m_file_name,*this,plines,72);

	//5.native_systemID
	put_Hollerith_string(m_native_systemID,*this,plines,72);

	//6.preprocessor_version
	put_Hollerith_string(m_preprocessor_version,*this,plines,72);
	
	//7.number_of_bits_of_integer
	put_integer(m_number_of_bits_of_integer,*this,plines,72);

	//8.magnitude_single_precision
	put_integer(m_magnitude_single_precision,*this,plines,72);

	//9.significance_single_precision
	put_integer(m_significance_single_precision,*this,plines,72);
	
	//10.magnitude_double_precision
	put_integer(m_magnitude_double_precision,*this,plines,72);

	//11.significance_double_precision
	put_integer(m_significance_double_precision,*this,plines,72);

	//12.productID_receiver
	put_Hollerith_string(m_productID_receiver,*this,plines,72);

	//13.model_space_scale
	put_real(m_model_space_scale,*this,plines,72);

	//14.unit_flag
	put_integer(m_unit_flag,*this,plines,72);

	//15.unit_name
	put_Hollerith_string(m_unit_name,*this,plines,72);

	//16.max_number_of_line_weight_gradations
	put_integer(m_max_number_of_line_weight_gradations,*this,plines,72);

	//17.width_of_max_line_weight
	put_real(m_width_of_max_line_weight,*this,plines,72);

	//18.DateTime_File_generation
	put_Hollerith_string(m_DateTime_File_generation,*this,plines,72);

	//19.min_resolution
	put_real(m_min_resolution,*this,plines,72);

	//20.max_coordinate_value
	put_real(m_max_coordinate_value,*this,plines,72);

	//21.author_name
	put_Hollerith_string(m_author_name,*this,plines,72);

	//22.author_organazation
	put_Hollerith_string(m_author_organazation,*this,plines,72);

	//23.version_flag
	put_integer(m_version_flag,*this,plines,72);

	//24.drafting_standard_flag
	put_integer(m_drafting_standard_flag,*this,plines,72);

	//25.DateTime_Model_generation
	put_Hollerith_string(m_DateTime_Model_generation,*this,plines,72);
	
	//26.application_protocolID
	put_Hollerith_string(m_application_protocolID,*this,plines,72);

	//Write out GSEC string to the file.
	int npl=plines.size();
	for(int j=0;j<npl;j++,lineID++){
		std::string& plj=*(plines[j]);
		ofsr<<std::left<<std::setw(72)<<std::setfill(' ')<<plj;
		ofsr<<'G'<<std::setw(7)<<std::setfill('0')<<std::right<<lineID<<std::endl;
	}
	
	return plines.size();	
}
