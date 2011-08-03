/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// mgSysGLList.cpp : mgSysGLList ÉNÉâÉXÇÃimplementationÅB
#include "MGCLStdAfx.h"
#include "mg/Curve.h"
#include "mgGL/SysGLList.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////Constructor/////////////

mgSysGLList::mgSysGLList(){;}

//Copy constructor.
mgSysGLList::mgSysGLList(
	const mgSysGLList& list2
):m_sysgls(list2.m_sysgls),m_sysgls_zebra(list2.m_sysgls_zebra){
	mgSysGLList* list2p=const_cast<mgSysGLList*>(&list2);
	list2p->clear();
}

//////////Destructor//////////////
mgSysGLList::~mgSysGLList(){clear();}

///////////////Operator overload/////////////////

//Assignment operator.
mgSysGLList& mgSysGLList::operator=(
	const mgSysGLList& list2
){
	clear();
	m_sysgls=list2.m_sysgls;
	m_sysgls_zebra=list2.m_sysgls_zebra;
	mgSysGLList* list2p=const_cast<mgSysGLList*>(&list2);
	list2p->clear();
	return *this;
}

///////////////Member fucntions/////////////////

//Clear the list.
void mgSysGLList::clear(){
	iterator i=m_sysgls.begin(), ie=m_sysgls.end();
	for(; i!=ie; i++){
		glDeleteLists((*i)->dlist_name(),1);
		//std::cout<<" clear delete Sysgl="<<(**i)<<std::endl;
		delete *i;
	}
	m_sysgls.clear();

	zebra_iterator j=m_sysgls_zebra.begin(), je=m_sysgls_zebra.end();
	for(; j!=je; j++){
		glDeleteLists((*j)->dlist_name(),1);
		//std::cout<<" clear delete Sysgl="<<(**j)<<std::endl;
		delete *j;
	}
	m_sysgls_zebra.clear();
}

//Delete all the display lists that have the fucntion_code fc.
//make_RC_current() of MGOpenGLView is necessary to use.
//function's retrurn value is true if any one is deleted.
bool mgSysGLList::delete_lists_by_function_code(size_t fc){
	bool deleted=false;
	iterator i=m_sysgls.begin(), in, ie=m_sysgls.end();
	in=i;
	while(i!=ie){
		in++;
		mgSysGL* sysgl=*i;
		if(sysgl->function_code()==fc){
			glDeleteLists(sysgl->dlist_name(),1);
			//std::cout<<" F delete Sysgl="<<(**i)<<std::endl;
			delete sysgl; deleted=true;
			m_sysgls.erase(i);
		}
		i=in;
	}

	zebra_iterator j=m_sysgls_zebra.begin(), jn, je=m_sysgls_zebra.end();
	jn=j;
	while(j!=je){
		jn++;
		mgSysZebraGL* sysgl=*j;
		if(sysgl->function_code()==fc){
			glDeleteLists(sysgl->dlist_name(),1);
			//std::cout<<" F delete Sysgl="<<(**j)<<std::endl;
			delete sysgl; deleted=true;
			m_sysgls_zebra.erase(j);
		}
		j=jn;
	}

	return deleted;
}

//Delete all the display lists that have the object_id oi.
//make_RC_current() of MGOpenGLView is necessary to use.
//function's retrurn value is true if any one is deleted.
bool mgSysGLList::delete_lists_by_object_id(
	const MGGel* oi,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be appended.
){
	bool deleted=false;
	iterator i=m_sysgls.begin(), in, ie=m_sysgls.end();
	in=i;
	while(i!=ie){
		in++;
		mgSysGL* sysgl=*i;
		if(sysgl->includes(oi)){
			glDeleteLists(sysgl->dlist_name(),1);
			functions.push_back(sysgl); 
			deleted=true;
			//std::cout<<" O delete Sysgl="<<(*sysgl)<<std::endl;
			m_sysgls.erase(i);
		}
		i=in;
	}

	zebra_iterator j=m_sysgls_zebra.begin(), jn, je=m_sysgls_zebra.end();
	jn=j;
	while(j!=je){
		jn++;
		mgSysGL* sysgl=*j;
		if(sysgl->includes(oi)){
			glDeleteLists(sysgl->dlist_name(),1);
			functions.push_back(sysgl); 
			deleted=true;
			//std::cout<<" O delete Sysgl="<<(*sysgl)<<std::endl;
			m_sysgls_zebra.erase(j);
		}
		j=jn;
	}

	return deleted;
}

//Delete the display list that have the fucntion_code fc and the object id gel.
//make_RC_current() of MGOpenGLView is necessary to use.
//function's retrurn value is true if any one is deleted.
bool mgSysGLList::delete_lists_by_function_object_code(
	size_t fc, const MGGel* gel
){
	bool deleted=false;

	iterator i=m_sysgls.begin(), in, ie=m_sysgls.end();
	in=i;
	while(i!=ie){
		in++;
		if((**i).function_code()==fc){
			if((**i).includes(gel)){
				glDeleteLists((*i)->dlist_name(),1);
				delete *i; deleted=true;
				//std::cout<<" O delete Sysgl="<<(**i)<<std::endl;
				m_sysgls.erase(i);
			}
		}
		i=in;
	}

	zebra_iterator j=m_sysgls_zebra.begin(), jn, je=m_sysgls_zebra.end();
	jn=j;
	while(j!=je){
		jn++;
		if((**j).function_code()==fc){
			if((**j).includes(gel)){
				glDeleteLists((*j)->dlist_name(),1);
				delete *j; deleted=true;
				//std::cout<<" O delete Sysgl="<<(**j)<<std::endl;
				m_sysgls_zebra.erase(j);
			}
		}
		j=jn;
	}

	return deleted;
}

//Draw all the objects by calling glCallList in this list.
//make_RC_current() of MGOpenGLView is necessary to use.
void mgSysGLList::draw_list()const{
	const_iterator i=m_sysgls.begin(), ie=m_sysgls.end();
	for(;i!=ie; i++)
		glCallList((*i)->dlist_name());
	const_zebra_iterator j=m_sysgls_zebra.begin(), je=m_sysgls_zebra.end();
	for(;j!=je; j++)
		glCallList((*j)->dlist_name());
}

//Test if this list includes the fucntion code fc's SysGL or not.
bool mgSysGLList::includes(size_t fc)const{
	const_iterator i=m_sysgls.begin(), ie=m_sysgls.end();
	for(;i!=ie; i++){
		if((**i).function_code()==fc)
			return true;
	}
	const_zebra_iterator j=m_sysgls_zebra.begin(), je=m_sysgls_zebra.end();
	for(;j!=je; j++){
		if((**j).function_code()==fc)
			return true;
	}

	return false;
}

//Invoke pre_transform_process of m_sysgls_zebra;
void mgSysGLList::pre_transform_process(){
	if(m_sysgls_zebra.empty())
		return;
	mgSysZebraGL* zebra=m_sysgls_zebra.back();
	zebra->pre_transform_process();
}

size_t mgSysGLList::push_back(size_t fc,const MGGel* oi){
	mgSysGL* sgl=new mgSysGL(fc,oi);
	m_sysgls.push_back(sgl);
	//std::cout<<"      Add Sysgl="<<(*sgl)<<std::endl;
	return sgl->dlist_name();
}
size_t mgSysGLList::push_front(size_t fc,const MGGel* oi){
	mgSysGL* sgl=new mgSysGL(fc,oi);
	m_sysgls.push_front(sgl);
	return sgl->dlist_name();
}

//sysgl must be a newed object, and the ownership will be 
//transfered to this.
size_t mgSysGLList::push_back(mgSysGL* sysgl){
	mgSysZebraGL* zebra=dynamic_cast<mgSysZebraGL*>(sysgl);
	if(zebra){
		m_sysgls_zebra.push_back(zebra);
	}else{
		m_sysgls.push_back(sysgl);
		//std::cout<<"      Add Sysgl="<<(*sysgl)<<std::endl;
	}
	return sysgl->dlist_name();
}
