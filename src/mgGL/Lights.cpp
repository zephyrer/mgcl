/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Light.h"
#include "mgGL/Lights.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////

//assignment
MGLights& MGLights::operator=(const MGLights& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	m_lights.clear();
	size_t n=gel2.m_lights.size();
	for(size_t i=0; i<n; i++){
		m_lights.push_back(gel2.m_lights[i]->clone());
	}
	return *this;
}
MGLights& MGLights::operator=(const MGGel& gel2){
	const MGLights* gel2_is_this=dynamic_cast<const MGLights*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGLights::operator<(const MGLights& gel2)const{
	return m_lights.size()<gel2.m_lights.size();
}
bool MGLights::operator<(const MGGel& gel2)const{
	const MGLights* gel2_is_this=dynamic_cast<const MGLights*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGLights* MGLights::clone()const{
	MGLights* lights=new MGLights;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++){
		lights->push_back(static_cast<MGLight*>((**i).clone()));
	}
	return lights;
}

void MGLights::exec()const{
	if(undefined()) return;

	size_t n=size();
	for(size_t i=0; i<n; i++) m_lights[i]->exec();
}

// Output function.
std::ostream& MGLights::out(std::ostream& ostrm) const{
	ostrm<<"Lgihts=";MGGLAttrib::out(ostrm);
	size_t n=size();
	ostrm<<",number of lights="<<n<<std::endl;
	for(size_t i=0; i<n; i++) m_lights[i]->out(ostrm);
	return ostrm;
}

//add one light.
//Function's return value is the numbe of lights defined.
size_t MGLights::push_back(MGLight* light){
	size_t lnum=m_lights.size();
	light->set_light_number(lnum);
	m_lights.push_back(light);
	return m_lights.size();
}

//Read all member data.
void MGLights::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	size_t n=size();
	buf>>n;
	for(size_t i=0; i<n; i++) push_back(static_cast<MGLight*>(buf.ReadPointer()));
}
//Write all member data
void MGLights::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	size_t n=size();
	buf<<n;
	const_iterator i=begin(), ie=end();	
	for(; i!=ie; i++) buf.WritePointer(*i);
}
