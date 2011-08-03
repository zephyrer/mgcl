/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// Copyright (c) 2011 Yuzi Mizuno and DG Technologies, Inc.
// 
// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without 
// limitation the rights to use, copy, modify, merge, publish, distribute, 
// sublicense, and/or sell copies of the Software, and to permit persons to 
// whom the Software is furnished to do so, subject to the following 
// conditions:
// 
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
	for(; first != last; ++first){
		if(first->second == factory){
			first = m_map.erase(first);
		}
	}
}
