/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __IGESOFSTREAM_H__)
#define __IGESOFSTREAM_H__

#include <fstream>
#include <string>
#include "mg/Pvector.h"
#include "mgIges/IgesFstream.h"
#include "mgIges/IgesDEStatusNumber.h"

// forward declarations
class MGObject;
class MGGroup;
class MGAttribedGel;
class MGGel;
class MGIgesPD;

/** @addtogroup FileInputOutput
 *  @{
 */

///MGOgesIfstream write out to *.iges file, transforming MGCL objects to IGES objects.
class MGIgesOfstream: public MGIgesFstream{

// Constructors.
public:

	/// Creates an object of class MGIgesOfstream with a filename.
	MGIgesOfstream(const char* filename=0);

	/// Destroys an object of class MGIgesOfstream.
	~MGIgesOfstream();

/// Operator overload.
public:

	/// Writes an object of the MGCL as an object of the IGES.
	MGIgesOfstream& operator<<(const MGGel&);

	/// Writes all objects of the MGGroup stored in group
	/// as IGES objects.
	MGIgesOfstream& operator<<(const MGGroup&);

/// Member functions.
public:

	///Get the reference of the output stream.
	std::ofstream& get_ofstream(){return m_ofstream;};

	///Test if output stream is open.
	bool is_open(){return m_ofstream.is_open();};

	///Open the file. When this is opened already, this will be closed, then will
	///be opened.
	void open(const char* filename);
	void close();

	///Create de as a newed object, and push back to this MGIgesOfstream.
	///The ownership of pd is transfered to the created DE, and the ownership of
	///the created DE will be transfered to this MGIgesOfstream.
	///Function's return value is the DE id created. If DE pointer is necessary,
	///get_DE() is the one.
	///SubordinateEntitySwitch is set as independent when gel is input, and
	///set as PLDependent when gel is not input(input as null).
	int create_de(
		MGIgesPD* pd,///<newed MGIgesPD object.
		const std::string& EntityLabel,
		int ses=0,///Subordinate entity switch. Default=independent.
		const MGAttribedGel* gel=0,
		int FormNumber=0	///<Form number
	);

	///Test if input stream is good().
	///True will be returned when good.
	bool good()const{return m_ofstream.good();};

	///Append a line of parameter data section record.
	///pline must be a newed object, and the ownership will be transfered to
	///this MGIgesOfstream class object.
	void append_param_line(MGIgesParamLine* pline){m_plines.push_back(pline);};

	///Get the next parameter data section record counter.
	int get_next_param_line_count()const{return m_plines.size()+1;};

	///get the output line number of Parameter Data Section.
	int get_line_number_of_PD()const{return m_plines.size();};

	///Write out start and terminate section into output file stream.
	void write_out_start_section();
	void write_out_terminate_section();
	
	///Write out DE lines to output file stream, together with PD lines.
	///write_out_DE_PD_lines outputs all the DEs and PDs using 
	///write_out_PD_pline for each DE.
	void write_out_DE_PD_lines();

	///Write out all the Paramete Data Lines in m_plines to output file stream.
	///m_plines[i] is one parameter data lines of IGES file
	///for 0<=i<m_plines.size().
	void write_out_PD_plines();

	///Obtain output file stream reference.
	std::ofstream& ofstrm(){return m_ofstream;};

private:

	std::ofstream m_ofstream;///<Output file stream.

	MGPvector<MGIgesParamLine> m_plines;
			///<vector of MGIgesParamLine, and controls
			///<line count of Paramete Data Section.
			///<m_plines[i-1] is line number i Parameter Data Section.

	///Initialize all the member data to the default value.
	void initialize(const char* filename=0);

	///Write out DEpointer's Paramete Data Lines into m_plines.
	///Paramete Data Lines of m_DrectoryEntries[DEpointer] will be
	///output into m_plines.
	///Function's return value is the number of lines output.
	int write_out_PD_pline(int DEpointer);

};

/** @} */ // end of FileInputOutput group
#endif // __IGESOFSTREAM_H__
