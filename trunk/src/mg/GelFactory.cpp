/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/GelFactory.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace
{
	bool g_valid = false; // �N���X�I�u�W�F�N�g���L�����ǂ���
}

// �R���X�g���N�^�[
MGGelFactoryRegistry::MGGelFactoryRegistry()
{
	g_valid = true;
}

// �f�X�g���N�^�[
MGGelFactoryRegistry::~MGGelFactoryRegistry()
{
	g_valid = false;
}

// �V���O���g���A�N�Z�X�B
MGGelFactoryRegistry* MGGelFactoryRegistry::get_instance()
{
	static MGGelFactoryRegistry registry;
	return &registry;
}

// �N���X�I�u�W�F�N�g���L�����ǂ�����Ԃ��B
bool MGGelFactoryRegistry::is_valid()
{
	return g_valid;
}

// ���O����R�}���h�I�u�W�F�N�g���쐬����B
MGGel* MGGelFactoryRegistry::create_gel(const KeyType& name) const
{
	TypeMap::const_iterator it = m_map.find(name);
	if(it == m_map.end()){
		//throw UnknownGelError(name);
		return 0;
	}

	return it->second->create_gel();
}

// �t�@�N�g���[�I�u�W�F�N�g��o�^����B
void MGGelFactoryRegistry::register_factory(
	const KeyType& name,
	MGGelFactoryBase* factory)
{
	m_map.insert(std::make_pair(name, factory));
}

// �t�@�N�g���[�I�u�W�F�N�g��o�^����폜����B
void MGGelFactoryRegistry::unregister_factory(const KeyType& name)
{
	m_map.erase(name);
}

// �t�@�N�g���[�I�u�W�F�N�g��o�^����폜����B
void MGGelFactoryRegistry::unregister_factory(MGGelFactoryBase* factory)
{
	TypeMap::iterator first = m_map.begin(), last = m_map.end();
	for(; first != last;){
		if(first->second == factory){
			first = m_map.erase(first);
		}
		else{
			++first;
		}
	}
}
