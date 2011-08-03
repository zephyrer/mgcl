/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "topo/Loop.h"
#include "topo/TrimLoop.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//MGTrimLoop express a loop and the start and end points are connected to
//which boundary of the loop.

//Copy constructor
MGTrimLoop::MGTrimLoop(const MGTrimLoop& linf2)
:m_loop(linf2.m_loop),m_start_loopid(linf2.m_start_loopid), m_start(linf2.m_start),
m_end_loopid(linf2.m_end_loopid),m_end(linf2.m_end){
	MGTrimLoop* linf2p=const_cast<MGTrimLoop*>(&linf2);
	linf2p->m_loop=0;
}

//Ordinary constructor.
MGTrimLoop::MGTrimLoop(
	MGLoop* loop,	//Newed loop pointer.
	int star_loop_id,MGLEPoint& start_lep, int end_loop_id,MGLEPoint& end_lep
):m_loop(loop),m_start_loopid(star_loop_id), m_start(start_lep),
		m_end_loopid(end_loop_id), m_end(end_lep){;}

//Destructor.
MGTrimLoop::~MGTrimLoop(){
	delete m_loop;
}

MGTrimLoop& MGTrimLoop::operator=(const MGTrimLoop& tloop2){
	delete m_loop;
	m_loop=tloop2.m_loop;
	m_start_loopid=tloop2.m_start_loopid;
	m_start=tloop2.m_start;

	m_end_loopid=tloop2.m_end_loopid;
	m_end=tloop2.m_end;
	MGTrimLoop* tloop2p=const_cast<MGTrimLoop*>(&tloop2);
	tloop2p->m_loop=0;
	return *this;
}

void MGTrimLoop::set_null(){
	delete m_loop;
	m_loop=0;
}

//Stream output of the content.
std::ostream& operator<< (std::ostream& ostrm, const MGTrimLoop& tloop){
	ostrm<<"MGTrimLoop="<<&tloop<<", m_loop="<<tloop.m_loop<<"="<<std::endl;
	ostrm<<*(tloop.m_loop)<<",";
	ostrm<<"   Start point="<<tloop.m_start_loopid;
	if(tloop.m_start_loopid<0 || tloop.m_start_loopid==4)
		ostrm<<", "<<tloop.m_start<<", "<<std::endl;
	else
		ostrm<<", "<<std::endl;

	ostrm<<"   End point="<<tloop.m_end_loopid;
	if(tloop.m_end_loopid<0 || tloop.m_end_loopid==4)
		ostrm<<", "<<tloop.m_end<<std::endl;
	else
		ostrm<<std::endl;

	return ostrm;
}

//Valid only when start_is_on_boundary();
size_t MGTrimLoop::start_loopid()const{
	size_t loopid=0;
	if(m_start_loopid<0)
		loopid=-m_start_loopid;

	return loopid;
}

size_t MGTrimLoop::end_loopid()const{
	size_t loopid=0;
	if(m_end_loopid<0)
		loopid=-m_end_loopid;

	return loopid;
}

MGLoop* MGTrimLoop::release_loop(){
	MGLoop* loop=m_loop;
	m_loop=0;
	return loop;
}

//Set used flag of used_loops@as true@for the both end loops id.
void MGTrimLoop::set_used_loop_flag(std::vector<bool>& used_loops)const{
	used_loops[end_loopid()]=true;
	used_loops[start_loopid()]=true;
}
