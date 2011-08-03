/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// MGOfstream.cpp: MGOfstream クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////

#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Gel.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////
//Find the input prt's map address.
//If found, position(pid) will be returned.
//If not found, null(0) will be returned.
//Here position means std::stremap of the m_file file where the ptr' is stored.
long MGOfstream::find(const MGGel* ptr){
	long pid = 0;
	mapitr mitr = m_map->pmap.find(ptr);
	if (mitr != m_map->pmap.end()){
		pid = (*mitr).second;
	}
	return pid;
}

//Insert the ptr into the map. Function's return value is:
//True: if ptr did not exist in the map and insertion succeeded.
//False: if ptr did exist in the map and insertion failed.
bool MGOfstream::insert(const MGGel* ptr, long pid){
	std::pair<mapitr, bool> ret = m_map->pmap.insert(pmap::value_type(ptr,pid));
	return ret.second;
}

void MGOfstream::mapClear(){m_map->pmap.clear();}

int MGOfstream::tellp(){
	return m_position;
}

// Pointer base のオブジェクトファイル出力関数
// 戻り値はオブジェクトのPID(識別ID)。
//This is an internal program. Ordinary users should not use this function.
//operator<< should be used instead.
long MGOfstream::WritePointer(const MGGel* gel){
	long pid;
	if(!gel){//When null pointer.
		pid=(*this).tellp();
		(*this)<<0x00000000L;
		return pid;
	}

	//自分自身が既にmapに登録されているかどうか調べる
	pid = (*this).find(gel);
	if(pid){
		//When gel is found in the map.
		(*this)<<pid;
		return pid;
	}

	(*this)<<0xffffffffL;//The header that indicates an object is followed.
	pid = (*this).tellp();	// PIDを決める(streamのPosition)
	(*this).insert(gel,pid);//mapに自分自身のPIDを登録

	long tid = gel->identify_type();
	(*this) << tid;
	gel->WriteMembers((*this));//データメンバの書き出し
	return pid;
}

//バイト列の並び順をUNIXスタイルに変更して出力する関数群
MGOfstream& MGOfstream::write1Byte(const void *ps){
	write(ps,1);
	return *this;
}

MGOfstream& MGOfstream::write2Byte(const void *ps2){
	char buf[2];
	char* tmp = (char*)ps2;
	buf[0]=tmp[1]; buf[1]=tmp[0];
	write(buf,2);
	return *this;
}

MGOfstream& MGOfstream::write4Byte(const void *ps4){
	char buf[4];
	char* tmp = (char*)ps4;
	buf[0]=tmp[3]; buf[1]=tmp[2]; buf[2]=tmp[1]; buf[3]=tmp[0];
	write(buf,4);
	return *this;
}

MGOfstream& MGOfstream::write8Byte(const void *ps8){
	char buf[8];
	char* tmp = (char*)ps8;
	buf[0]=tmp[7]; buf[1]=tmp[6]; buf[2]=tmp[5]; buf[3]=tmp[4];
	buf[4]=tmp[3]; buf[5]=tmp[2]; buf[6]=tmp[1]; buf[7]=tmp[0];
	write(buf,8);
	return *this;
}

MGOfstream& MGOfstream::writenChar(char* ps, size_t n){
	write(ps,n);
	return *this;
}
