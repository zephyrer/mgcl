/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Attrib.h"
#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/Group.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGGel
// Implementation of MGGel.
//MGGel is an abstract class which represents a group element.
//Gel is the abbreviation of group element.
//Subclasses of MGGel are:
//(1) MGAttribedGel, or (2) MGAttrib.
//MGGel provides functions of serialization of objects.
//All the objects of MGGel subclasses can be serialized using
//MGGroup::make_file(), and MGGroup constructor.

//Virtual Destructor
MGGel::~MGGel(){;}

//Construct a null newed MGGel from the type id TID.
//Object handled by MGIfstream or MGOfstream is only the following objects.
MGGel* MGNullGel(long TID){
	long id=TID; id&=0xff000000;
	switch(id){
		case MGOBJECT_TID:	return MGNullObj(TID);
		case MGGROUP_TID:	return MGNullGroup(TID);
		case MGATTRIB_TID:	return MGNullAttrib(TID);
		//case MGATTRIBEDOBJ_TID:return MGNullAttribedObject(TID);
	}
	return 0;
}

bool MGGel::operator<(const MGGel& gel2)const{
	long id1=identify_type();
	long id2=gel2.identify_type();
	return id1<id2;
}

//Read all member data.
void MGGel::ReadMembers(MGIfstream& buf){;}
//Write all member data
void MGGel::WriteMembers(MGOfstream& buf)const{;}

//Output the content as std::string.
//The output string is the same as cout<<MGGel.
std::string MGGel::string_content()const{
	std::ostringstream s;
#ifdef WIN32
	//added by Tomoko.
	//Global Localeをセットすること。
	s.imbue(std::locale::empty());
#endif
	out(s);
	return s.str();
}

//Determine if this is one of the input types or not.
//Function's return value is true if this is one of the input types.
bool MGGel::type_is(const MGAbstractGels& types)const{
	MGAbstractGels::const_iterator i=types.begin(), ie=types.end();
	for(; i!=ie; i++){
		MGGEL_KIND agel=(*i).first;
		long tid1=identify_type()&agel;
		long tid2=(*i).second&agel;
		if(tid1==tid2)
			return true;
	}
	return false;
}

//////////// MGGel ////////////
ostream& operator<< (ostream& ostrm, const MGGel& gel){
	return gel.out(ostrm);
}
