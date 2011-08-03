/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGEDESTATUSNUMBER_H__)
#define __MGIGEDESTATUSNUMBER_H__

#include <string>

///MGIgesDEStatusNumber describes the Status Number of a directory entry section.
class MGIgesDEStatusNumber{
// Constructors.
public:

	///Subordinate entity switch
	enum SESwitch{
		independent=0,
		PDependent=1,	///<Physically dependent
		LDependent=2,	///<Logically dependent
		PLDependent=3	///<Both Physically and Logically dependent
	};

	/// Constructs an object of class MGIgesDEStatusNumber.
	///Default constructor, includes all the defalut value of MGCL.
	MGIgesDEStatusNumber();

	MGIgesDEStatusNumber(
		short BlankStatus,///<<0:Visible, 1:Blanked;
		SESwitch SubordinateEntitySwitch,
				///<<0:Independent, 1:Physically Dependent, 2: Logically Dependent, 3:Both 1 and 2.
		short EntityUseFlag,
				///<0:Geometry, 1:Annotation, 2:Definition, 3:Other,
				///<4:Logical/Positional, 5:2D Parametric, 6:Construction Geometry;
		short Hierarchy///<0:Global top down, 1:Global defer, 2:Use hierarchy property;
	);

	short blankStatus()const{return m_BlankStatus;};
	SESwitch subordinateEntitySwitch()const{return (SESwitch)m_SubordinateEntitySwitch;};
	short entityUseFlag()const{return m_EntityUseFlag;};
	short hierarchy()const{return m_Hierarchy;};
	void read_in(const std::string& status);
	void set_as_blank(){m_BlankStatus=1;};
	void set_SubordinateEntitySwitch(SESwitch eswitch){m_SubordinateEntitySwitch=eswitch;};

private:
	short m_BlankStatus;///<0:Visible, 1:Blanked;
	short m_SubordinateEntitySwitch;
				///<0:Independent, 1:Physically Dependent, 2: Logically Dependent, 3:Both 1 and 2.
	short m_EntityUseFlag;
				///<0:Geometry, 1:Annotation, 2:Definition, 3:Other,
				///<4:Logical/Positional, 5:2D Parametric, 6:Construction Geometry;
	short m_Hierarchy;///<0:Global top down, 1:Global defer, 2:Use hierarchy property;
};

#endif // __MGIGEDESTATUSNUMBER_H__
