/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined(MGINCLUDEGUARD_GELFACTORY__)
#define MGINCLUDEGUARD_GELFACTORY__
#include <map>
#include "mg/MGCL.h"

class MGGel;

/// @class MGGelFactoryBase
/// Factory Method �p�̃C���^�[�t�F�C�X
struct MGGelFactoryBase
{
	/// �T�u�N���X�Ŏ������� MGGel �N���X�̃I�u�W�F�N�g��Ԃ��B
	/// @return MGGel new �I�u�W�F�N�g�B
	virtual MGGel* create_gel() const = 0;
};

/// @class MGGelFactoryT
/// @sa MGGelFactoryBase
template <typename T>
struct MGGelFactoryT : public MGGelFactoryBase
{
	/// @return MGGel new �I�u�W�F�N�g�B
	virtual T* create_gel() const
	{
		return new T;
	}
};

/// @class MGGelFactoryRegistry
/// �I�u�W�F�N�g�t�@�N�g���[�N���X�B
///
/// �T���v��
///
/// MGObject* MGNullObj(long TID){
///     MGGelFactoryRegistry* reg = MGGelFactoryRegistry::get_instance();
///     return static_cast<MGObject*>(reg->create_gel(TID));
/// }
class MGCLASS MGGelFactoryRegistry
{
public:
	typedef long KeyType;
	typedef std::map<KeyType, MGGelFactoryBase*> TypeMap;

	/// �f�X�g���N�^�[
	/// @post is_valid() �� false ��Ԃ��B
	/// @throws n/a
	~MGGelFactoryRegistry();

	/// �V���O���g���A�N�Z�X�B
	/// @return MGGelFactory �t�@�N�g���[�I�u�W�F�N�g�B
	/// @post is_valid() �� true ��Ԃ��B
	/// @throws n/a
	static MGGelFactoryRegistry* get_instance();

	/// �N���X�I�u�W�F�N�g���L�����ǂ�����Ԃ��B
	/// @return bool �L���Ȃ�΁i�f�X�g���N�^�[���Ă΂�Ă��Ȃ���΁jtrue
	/// @throws n/a
	static bool is_valid();

	/// ���O����I�u�W�F�N�g���쐬����B
	/// @param[in] name ���O
	/// @return MGGel new �I�u�W�F�N�g�B
	/// @throws n/a
	/// ���݂��Ȃ��^�C�v�̏ꍇ�A�k����Ԃ��B
	MGGel* create_gel(const KeyType& name) const;

	/// �t�@�N�g���[�I�u�W�F�N�g��o�^����B
	/// @param[in] name ���O
	/// @param[in] factory �t�@�N�g���[�I�u�W�F�N�g
	/// @pre is_valid() �� true ��Ԃ��B
	void register_factory(
		const KeyType& name,
		MGGelFactoryBase* factory);

	/// �t�@�N�g���[�I�u�W�F�N�g��o�^����폜����B
	/// @param[in] name ���O
	/// @pre is_valid() �� true ��Ԃ��B
	/// @throws n/a
	void unregister_factory(const KeyType& name);

	/// �t�@�N�g���[�I�u�W�F�N�g��o�^����폜����B
	/// @param[in] factory �t�@�N�g���[�I�u�W�F�N�g
	/// @pre is_valid() �� true ��Ԃ��B
	/// @throws n/a
	void unregister_factory(MGGelFactoryBase* factory);

private:
	/// �R���X�g���N�^�[
	MGGelFactoryRegistry();

	MGGelFactoryRegistry(const MGGelFactoryRegistry& other);
	MGGelFactoryRegistry& operator=(const MGGelFactoryRegistry& other);

	TypeMap m_map; ///< �I�u�W�F�N�g���ƃI�u�W�F�N�g�t�@�N�g���[�̎����B
};

/// @class MGAutoGelRegister
/// �t�@�N�g���[���W�X�g���[�ɃG���g���[����֗��N���X
/// @sa AUTO_GEL_REGISTER
template <typename T>
class MGAutoGelRegister
{
	MGGelFactoryRegistry* m_pRegistry; /// Singleton �I�u�W�F�N�g�ւ̃|�C���^�[
	MGGelFactoryT<T> m_factory; /// T �I�u�W�F�N�g�̂��߂̃t�@�N�g���[

public:
	/// �R���X�g���N�^�[
	/// @param[in] name �I�u�W�F�N�g�̖��O
	MGAutoGelRegister(const MGGelFactoryRegistry::KeyType& name)
		 : m_pRegistry(MGGelFactoryRegistry::get_instance())
	{
		m_pRegistry->register_factory(name, &m_factory);
	}

	/// �f�X�g���N�^�[
	~MGAutoGelRegister()
	{
		if(MGGelFactoryRegistry::is_valid()){
			m_pRegistry->unregister_factory(&m_factory);
		}
	}
};

/// For internal use.
#define GEL_REG_JOIN( s1, s2 ) GEL_REG_JOIN__( s1, s2 )
/// For internal use.
#define GEL_REG_JOIN__( s1, s2 ) s1##s2
/// For internal use.
#define GEL_REG_MAKE_UNIQUE_NAME( prefix ) GEL_REG_JOIN( prefix, __LINE__ )

/// �V�K�� MG �̃T�u�N���X���`���A������V���A���C�Y�ΏۂƂ������ꍇ�́A
/// ���̎菇�𓥂ށB
///
/// 1. virtual long identify_type() const ���I�[�o�[���C�h���A
///    �K�؂Ȓl��Ԃ��悤��������B
///    mg/types.h �� enum MGGEL_TID �̋K���ɏ]�����ƁB
///
/// 2. �T�u�N���X�̎����t�@�C�� (cpp) �ɁA���̊֐��^�}�N���Ăяo����
///    ��s�L�q���邱�ƁB
///
///    AUTO_GEL_REGISTER(MGMyClass, 0x????????L);
///
///    �����ŁA0x????????L �� MGMyClass::identify_type �ƈ�v����l�Ƃ��邱�ƁB
#define AUTO_GEL_REGISTER(classname, gelname) \
	static MGAutoGelRegister<classname> \
		GEL_REG_MAKE_UNIQUE_NAME( glreg__ ) (gelname);

#endif // defined(MGINCLUDEGUARD_GELFACTORY__)
