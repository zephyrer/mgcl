/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//////////////////////////////////////////////////////////////////////
// MGIfstream.cpp: MGIfstream クラスのインプリメンテーション
//////////////////////////////////////////////////////////////////////

#include "MGCLStdAfx.h"
#include "mg/Ifstream.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//Find the input pid's map address.
//If found, MGGel* will be returned.
//If not found, null pointer(0) will be returned.
MGGel* MGIfstream::find(long pid){
	MGGel* ptr = NULL;
	mapitr mitr = m_map->pmap.find(pid);
	if (mitr != m_map->pmap.end()){
		ptr = (*mitr).second;
	}
	return ptr;
}

//Insert the ptr into the map. Function's return value is:
//True: if ptr did not exist in the map and insertion succeeded.
//False: if ptr did exist in the map and insertion failed.
bool MGIfstream::insert(long pid,MGGel* ptr){
	std::pair<mapitr, bool> ret = m_map->pmap.insert(pmap::value_type(pid,ptr));
	return ret.second;
}

//Clear map
void MGIfstream::mapClear(){m_map->pmap.clear();}

//バイト列の並び順をUNIXスタイルに変更して出力する関数群
MGIfstream& MGIfstream::read1Byte(void *ps){
	read(ps,1);
	return *this;
}

//バイト列の並び順をUNIXスタイルに変更して出力する関数群
MGIfstream& MGIfstream::read2Byte(void *ps2){
	char buf[2];
	read(buf,2);
	char* tmp = (char*)ps2;
	tmp[0]=buf[1]; tmp[1]=buf[0];
	return *this;
}

MGIfstream& MGIfstream::read4Byte(void *ps4){
	char buf[4];
	read(buf,4);
	char* tmp = (char*)ps4;
	tmp[0]=buf[3]; tmp[1]=buf[2]; tmp[2]=buf[1]; tmp[3]=buf[0];
	return *this;
}

MGIfstream& MGIfstream::read8Byte(void *ps8){
	char buf[8];
	read(buf,8);
	char* tmp = (char*)ps8;
	tmp[0]=buf[7]; tmp[1]=buf[6]; tmp[2]=buf[5]; tmp[3]=buf[4];
	tmp[4]=buf[3]; tmp[5]=buf[2]; tmp[6]=buf[1]; tmp[7]=buf[0];
	return *this;
}

MGIfstream& MGIfstream::readnChar(char* ps, size_t n){
	read(ps,n);
	return *this;
}

long MGIfstream::tellg(){
	return m_position;
}

