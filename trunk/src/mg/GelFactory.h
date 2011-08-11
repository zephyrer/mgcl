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
/// Factory Method 用のインターフェイス
struct MGGelFactoryBase
{
	/// サブクラスで実装して MGGel クラスのオブジェクトを返す。
	/// @return MGGel new オブジェクト。
	virtual MGGel* create_gel() const = 0;
};

/// @class MGGelFactoryT
/// @sa MGGelFactoryBase
template <typename T>
struct MGGelFactoryT : public MGGelFactoryBase
{
	/// @return MGGel new オブジェクト。
	virtual T* create_gel() const
	{
		return new T;
	}
};

/// @class MGGelFactoryRegistry
/// オブジェクトファクトリークラス。
///
/// サンプル
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

	/// デストラクター
	/// @post is_valid() が false を返す。
	/// @throws n/a
	~MGGelFactoryRegistry();

	/// シングルトンアクセス。
	/// @return MGGelFactory ファクトリーオブジェクト。
	/// @post is_valid() が true を返す。
	/// @throws n/a
	static MGGelFactoryRegistry* get_instance();

	/// クラスオブジェクトが有効かどうかを返す。
	/// @return bool 有効ならば（デストラクターが呼ばれていなければ）true
	/// @throws n/a
	static bool is_valid();

	/// 名前からオブジェクトを作成する。
	/// @param[in] name 名前
	/// @return MGGel new オブジェクト。
	/// @throws n/a
	/// 存在しないタイプの場合、ヌルを返す。
	MGGel* create_gel(const KeyType& name) const;

	/// ファクトリーオブジェクトを登録する。
	/// @param[in] name 名前
	/// @param[in] factory ファクトリーオブジェクト
	/// @pre is_valid() が true を返す。
	void register_factory(
		const KeyType& name,
		MGGelFactoryBase* factory);

	/// ファクトリーオブジェクトを登録から削除する。
	/// @param[in] name 名前
	/// @pre is_valid() が true を返す。
	/// @throws n/a
	void unregister_factory(const KeyType& name);

	/// ファクトリーオブジェクトを登録から削除する。
	/// @param[in] factory ファクトリーオブジェクト
	/// @pre is_valid() が true を返す。
	/// @throws n/a
	void unregister_factory(MGGelFactoryBase* factory);

private:
	/// コンストラクター
	MGGelFactoryRegistry();

	MGGelFactoryRegistry(const MGGelFactoryRegistry& other);
	MGGelFactoryRegistry& operator=(const MGGelFactoryRegistry& other);

	TypeMap m_map; ///< オブジェクト名とオブジェクトファクトリーの辞書。
};

/// @class MGAutoGelRegister
/// ファクトリーレジストリーにエントリーする便利クラス
/// @sa AUTO_GEL_REGISTER
template <typename T>
class MGAutoGelRegister
{
	MGGelFactoryRegistry* m_pRegistry; /// Singleton オブジェクトへのポインター
	MGGelFactoryT<T> m_factory; /// T オブジェクトのためのファクトリー

public:
	/// コンストラクター
	/// @param[in] name オブジェクトの名前
	MGAutoGelRegister(const MGGelFactoryRegistry::KeyType& name)
		 : m_pRegistry(MGGelFactoryRegistry::get_instance())
	{
		m_pRegistry->register_factory(name, &m_factory);
	}

	/// デストラクター
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

/// 新規に MG のサブクラスを定義し、それをシリアライズ対象としたい場合は、
/// 次の手順を踏む。
///
/// 1. virtual long identify_type() const をオーバーライドし、
///    適切な値を返すよう実装する。
///    mg/types.h の enum MGGEL_TID の規則に従うこと。
///
/// 2. サブクラスの実装ファイル (cpp) に、次の関数型マクロ呼び出しを
///    一行記述すること。
///
///    AUTO_GEL_REGISTER(MGMyClass, 0x????????L);
///
///    ここで、0x????????L は MGMyClass::identify_type と一致する値とすること。
#define AUTO_GEL_REGISTER(classname, gelname) \
	static MGAutoGelRegister<classname> \
		GEL_REG_MAKE_UNIQUE_NAME( glreg__ ) (gelname);

#endif // defined(MGINCLUDEGUARD_GELFACTORY__)
