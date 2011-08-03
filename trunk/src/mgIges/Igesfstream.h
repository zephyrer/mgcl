/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __IGESFSTREAM_H__)
#define __IGESFSTREAM_H__

#include "mg/Pvector.h"
#include "mgIges/Iges.h"
#include "mgIges/IgesGSec.h"
#include "mgIges/IgesDirectoryEntry.h"
#include "mgIges/IgesParamLine.h"

/** @addtogroup FileInputOutput
 *  @{
 */

///MGIgesFstream is a super class for MGIfstream and MGOfstream.
///MGIgesFstream holds the data for IGES, and provides a common functions
///to write and read IGES data.
class MGIgesFstream{

/// Constructors.
public:

	///Destructor;
	virtual ~MGIgesFstream(){;}

	///Initialize all the member data to the state of no_value_holding.
	virtual void initialize(const char* filename=0);

///Function's return value is the directory entry pointer pushed back.
	int push_back_DE(MGIgesDirectoryEntry* de);

	///Return directory entry point of DEid.
	MGIgesDirectoryEntry* directoryEntry(int DEid){
		return m_DirectoryEntries[DEid];
	};
	const MGIgesDirectoryEntry* directoryEntry(int DEid)const{
		return m_DirectoryEntries[DEid];
	};

	void clearStartSection(){m_StartSection=std::string();};
	void clearGSection(){m_GSection=MGIgesGSec();};
	void clearDirectoryEntries(){m_DirectoryEntries.clear();};
	void clear(){clearStartSection(); clearDirectoryEntries(); clearDirectoryEntries();};
	void set_initial_StartSection();

	const MGIgesGSec& GSection()const{return m_GSection;};
	MGIgesGSec& GSection(){return m_GSection;};

	///get the output line number of Start Section.
	int get_line_number_of_SS()const{return m_StartSection.size()/72+1;};

	///get the output line number of Global Sections.
	int get_line_number_of_GS()const{return m_nlineGSec;};

	///get the output line number of Directory Entries.
	int get_line_number_of_DE()const{return (m_DirectoryEntries.size()-1)*2;};

protected:

// Data members.

	std::string m_StartSection;///<Start section string data.
	MGIgesGSec m_GSection;///<Global section data.
	int m_nlineGSec;	///<Number of the lines of Global section.
	MGPvector<MGIgesDirectoryEntry> m_DirectoryEntries;///<Directry entry data vector.
				///<One pair of directory entry lines are stored in m_DirectryEntry[i].
				///<m_DirectoryEntries[0] is a dummy entry and has no meaning.
};

/** @} */ // end of FileInputOutput group
#endif // __IGESFSTREAM_H__
