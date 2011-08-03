/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//////////////////////////////////////////////////////////////////////
// MGOfstream.cpp: MGOfstream クラスのインプリメンテーション
//////////////////////////////////////////////////////////////////////
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Object.h"
#include "mg/Group.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using std::ofstream;
using std::ios;

//////////////////////////////////////////////////////////////////////
// Constructor and destructor.
//////////////////////////////////////////////////////////////////////

MGOfstream::MGOfstream()
:m_file(0), m_map(new MGOutPtrMap()),m_position(0), m_map_clear(true){}

MGOfstream::MGOfstream(const char *file, bool map_clear)
:m_file(0),m_map(new MGOutPtrMap()),m_position(0){
	open(file,map_clear);
}

MGOfstream::~MGOfstream(){
	close();
	delete m_map;
}

///////////////Operator overload//////////////////

//Write out an object to ofs. Written objects by this fucction are able
//to read by one of the following global function:
// 1. MGIfstream& operator>>(MGIfstream& file, MGGel*& obj);
// 2. operator>>(MGIfstream& file, MGPvector<MGGel>& vecObj);
// 3. operator>>(MGIfstream& file, MGPlist<MGGel>& listObj);
MGOfstream& MGOfstream::operator<< (const MGGel& obj){
	if(m_map_clear) mapClear();
	(*this)<<0xffffffffL;//Flag that indicates an object is followed.
	long pid = tellp();  // PIDを決める(streamのPosition)
	insert(&obj,pid);	 //mapに自分自身のPIDを登録

	//データメンバの書き出し
	long tid = obj.identify_type();
	(*this)<<tid;
	obj.WriteMembers(*this);
	return (*this);
}

/*void MGOfstream::operator<<(MGPvector<MGGel>& vecObj){
	std::vector<MGGel*>::const_iterator i=vecObj.begin(), ie=vecObj.end();
	for(; i!=ie; i++) (*this)<<(**i);
}*/

void MGOfstream::operator<<(const MGGroup& group){
	MGGroup::const_iterator i=group.begin(), ie=group.end();
	for(; i!=ie; i++) (*this)<<(**i);
}

/*
//The following two functions are provided to hold the compatibility with
//version 5.
void MGOfstream::operator<<(MGPvector<MGObject>& vecgel){
	MGPvector<MGObject>::const_iterator i=vecgel.begin(), ie=vecgel.end();
	for(; i!=ie; i++) (*this)<<(**i);
}
void MGOfstream::operator<<(MGPlist<MGObject>& listgel){
	MGPlist<MGObject>::const_iterator i=listgel.begin(), ie=listgel.end();
	for(; i!=ie; i++) (*this)<<(**i);
}*/

///////////////Member function//////////////////

void MGOfstream::close(){
	if(m_file) delete m_file;
	m_file=0;
}

//Function's return value is:
//=0: open succeeded.
//=1: file not found, or could not be opened.
//=2: file found, but, the format is not MGCL format.
int MGOfstream::open(const char* file, bool map_clear){
	m_map_clear=map_clear;
	int error=0;
	std::ofstream* file_new=0;
	file_new = new ofstream(file, ios::out | ios::binary);

	if(file_new && file_new->is_open()){
		if(m_file){
			m_file->close();
			delete m_file;
		}
		m_file=file_new;
		write(MGCL_Version(),8);
		write(MGCL_File_validity(),24);
		error=0;
	}else{
		delete file_new; file_new=0;
		error=1;
	}
	return error;
}

//Write out n bytes date in th buffer ps.
//This read data is row data, and the sequence will not be changed
//like read nByte.
void MGOfstream::write(const void* ps, size_t n){
	m_file->write((char*)ps,n);
	m_position+=n;
}
