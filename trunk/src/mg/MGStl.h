/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGStl_HH_
#define _MGStl_HH_

#include <map>
#include "mg/object.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/FSurface.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "Tl/TLData.h"

class mgTLData;

/** @addtogroup MGObjectRelated
 *  @{
 */

//Define MGStl Class.

///MGStl is a concrete class which represents an STL information.
///このクラスは複数の三角形を配列として保持する.
///全ての三角形の頂点の座標はファイルから読み込んだ順で、重複を取り除いて配列に格納される.
///m_indicesにはファイルから読み込んだ順番で各三角形の頂点の並びが格納されており
///その各要素には該当する頂点座標の配列の添え字が格納される.
///例えば、i番目の三角形の各頂点の座標の添え字はindices[i*3]、indices[i*3+1]、indices[i*3+2]
///という以上の順序で格納されている.
class MGCLASS MGStl : public MGObject{

	/// mapには1つのキーに対して複数の値は収められないため
	/// 三角形を構成する頂点のIDを収めた構造体を利用する
	class vertId{
	public:
		int id1;
		int id2;
		int id3;
	};

	/// 関数オブジェクトを格納する構造体
	class positionComp{
	public:
		bool operator()(const MGPosition& p1, const MGPosition& p2) const{
			/// lexicographical_compare style
			if(p1(0) < p2(0)){
				return true;
			}else if(p2(0) < p1(0)){
				return false;
			}
			if(p1(1) < p2(1)){
				return true;
			}else if(p2(1) < p1(1)){
				return false;
			}
			if(p1(2) < p2(2)){
				return true;
			}else if(p2(2) < p1(2)){
				return false;
			}
			/// p1 == p2
			return false;
		}		
	};

	/// 三角形の頂点の座標、頂点の番号を格納するmapの別名を定義.
	typedef std::map<MGPosition, int, positionComp> triangleMap;

	/// 頂点のインデックス、頂点の座標を格納するmapの別名を定義.
	typedef std::map<int, MGPosition> IndexPosMap;

	/// 三角形のインデックス、3つの頂点のインデックスを格納するmapの別名を定義.
	typedef std::map<int, vertId> TriVertMap;

public:
	/// デフォルトコンストラクタ
	MGStl(void){;};

	/// コピーコンストラクタ
	MGStl(const MGStl& stl);

	///conversion constructor from tessellation data.
	MGStl(const mgTLDataVector& tlDataVector);
	MGStl(const mgTLData& tlData);

	///Constructor from triangle data, index+vertices.
	///This constructor uses the wc_zero to identify
	///different two positions as the same input position.
	MGStl(
		int nTriang, /// 三角形の数
		const int* triang, /// trang[i], [i+1], [i+2] make a triangle for i=0,...,nTriang-1
		const double* verts /// 頂点の座標
	);

	/// 仮想デストラクタ.
	~MGStl(void);

	/// 全ての頂点とボックスの座標に指定したMGVectorの値を加算.
	MGStl& operator+=(const MGVector& v);

	/// 全ての頂点とボックスの座標から指定したMGVectorの値を引く.
	MGStl& operator-=(const MGVector& v);

	/// 全ての頂点とボックスの座標に指定した値を掛ける.
	MGStl& operator*=(double scale);

	/// 与えられた変換を行う.
	MGStl& operator*=(const MGMatrix& mat);

	/// 与えられた変換によるトランスフォームを行う.
	MGStl& operator*=(const MGTransf& tr);

	/// 2つのMGStlオブジェクトが等しいか判定する.
	bool operator==(const MGStl& stl);

	/// 2つのMGStlオブジェクトを代入する.
	MGStl& operator=(const MGStl& stl);

	std::ostream& out(std::ostream& ostrm)const;

	/// STLファイルの図形のbounding boxを取得する.
	/// 戻り値：STLファイルの図形のbounding box
	const MGBox& box()const{return m_box;};

	///Construct new object by copying to newed area.
	///User must delete this copied object by "delete".
	MGStl* clone()const;

	/// Return This object's typeID.
	long identify_type()const;

	///Draw 3D curve in world coordinates.
	///The object is converted to curve(s) and is drawn.
	void drawWire(
		double span_length,	///<Line segment span length.
		int line_density=1	///<line density to draw a surface in wire mode.
	)const;

	///Draw 3D point(vertex) in world coordinates.
	///The object is converted to point(s) and is drawn.
	///This is valid only for topology objects or MGPoint.
	void draw3DVertex()const;

	///Shade the object in world coordinates.
	void shade(
		double span_length	///<Line segment span length.
	)const;

	///Compute the intersections of two objects.
	///***********MGStl does not support the intersection functions.******
	MGisects intersection(const MGObject& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGCurve& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGFSurface& obj2)const{return MGisects(this,obj2.object_pointer());};
	MGisects intersection(const MGSurface& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGFace& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGShell& obj2)const{return MGisects(this,&obj2);};

	///Get manifold dimension.
	unsigned manifold_dimension()const{return 2;};

	///Write all member data.
	void WriteMembers(MGOfstream& buf)const;

	///Read all member data.
	void ReadMembers(MGIfstream& buf);

	///return the ref to m_vecPos.
	const std::vector<MGPosition>& positions()const{return m_vecPos;};
	std::vector<MGPosition>& positions(){return m_vecPos;};

	///return the ref to m_vecPos.
	const std::vector<MGUnit_vector>& normals()const{return m_vecNormlTriang;};
	std::vector<MGUnit_vector>& normals(){return m_vecNormlTriang;};

	/// 三角形の個数を取得する.
	/// 戻り値：三角形の個数.
	int GetTriangleCount()const{return m_vecNormlTriang.size();};

	/// 指定した三角形の頂点座標が入れられている配列のインデックスを取得する.
	void GetVertIndices(
		int i, ///< [in]：三角形のインデックス(0 <= i < GetTriangleCount)
		size_t pos[3] ///< [out]：指定した三角形の各頂点座標の配列の添え字[i, j, k]
	)const;

	/// 引数で指定したパスのSTLファイルを読み込みメンバに値を設定する.
	/// またMGTolerance::wc_zeroに値を設定する.
	/// 戻り値: =0ファイルの読み込が成功
	///		  !=0 読み込まなかった。または失敗した(std::ifstreamのエラーコード）
	/// 事後条件:m_vecPos, m_vecNorml, m_indices, m_boxに図形の情報が格納される.
	///			MGTolerance::wc_zeroに値が設定される.
	int LoadFile(
		const char* strFilePath ///< [in]:読み込むSTLファイルへのパス
	);

	/// 指定されたパスにAscii形式のSTLファイルを保存する.
	/// 戻り値: =0ファイルの書き込みが成功
	///		  !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
	/// 事後条件:rSTLFilePathに指定したパスにAscii形式のSTLファイルが保存される.
	int SaveAscii(
		std::ofstream& fout
	)const;

	/// 指定されたパスにBinary形式のSTLファイルを保存する.
	/// 戻り値: =0ファイルの書き込みが成功
	///		  !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
	/// 事後条件:rSTLFilePathに指定したパスにBinary形式のSTLファイルが保存される.
	///仕様変更：引数にストリームを渡されるように変更09/09/17.
	int SaveBinary(
		const char* rSTLFilePath ///< [in]:ファイルの保存先のパス
	)const;

	/// 3点から作成した三角形の情報をメンバ変数に値を設定する.
	/// 事前条件:入力するpos1, pos2, pos3は反時計回りになっていること.
	/// 事後条件:m_vecPos, m_vecNormlTriang, m_indicesに三角形の情報が設定される.
	void push_back_triangle(
		const MGPosition& pos1,
		const MGPosition& pos2,
		const MGPosition& pos3,
		triangleMap& VertexMap
	);

	/// 三角形ごとに法線を表示する.
	void display_arrows()const;

	/// メンバーデータを直接セットする.
	void set_all_data(
		const std::vector<MGPosition>& vertices,
		const std::vector<int>& indices);

	/// Update the bounds.
	void update_bounds();

	/// Update the all normals of the triangles.
	void update_normals();

	/// オブジェクト名を返却する.
	std::string whoami()const{return "Stl";};

private:
	/// 一度ファイルを読み、メンバに値が設定された状態で再びファイルを読むと.
	/// メンバの値が上書きされる。これを避けるため、メンバの値を全てリセットする.
	void Initialize();

	/// 座標値をm_vecPosへpush_backし、m_box値をexpand()する.
	/// 事後条件:m_vecPosに要素が1つ追加される.
	void push_back_Position(
		const MGPosition& pos ///< [in]:m_vecPosに追加する頂点の座標.
		);

	/// 三角形を構成する頂点IDをm_indicesにpush_backし
	/// 指定した三角形について面法線ベクトルを計算しm_vecNormlTriangへpush_backする.
	/// 事前条件:m_vecPosに値が格納されている.
	/// 事後条件:１つの三角形の情報がm_vecNormlTriang、m_indecesに格納される.
	void add_indices_and_calc_normal(
	 const int vertId[3], ///< [in]:１つの三角形の頂点IDの配列.
	 int triIndex ///< [in]:法線を求める対象の三角形のインデックス
				  ///< (0 <= triIndex < GetTriangleCount).
	 );

	/// Ascii形式のファイルを読み込み全ての座標値を取得する.
	/// また図形のボックス枠を設定する.
	/// 戻り値: =0ファイルの読み込が成功
	///		  !=0 読み込みに失敗した(std::ifstreamのエラーコード)
	/// 事前条件:既にオープンされたファイルストリームを渡す.
	/// 事後条件:vecPosに全ての座標値が格納され、 boxにボックス枠の座標値が設定される.
	/// ファイルストリームが進む.
	int LoadAscii(
		std::ifstream& in, ///< [in/out]:既にオープンされたファイルストリーム .
		std::vector<MGPosition>& vecPos ///<[out]:ファイルから読み込んだ座標値の配列.
		);

	/// Binary形式のファイルを読み込み全ての座標値を取得する.
	/// また図形のボックス枠を設定する.
	/// 戻り値: =0ファイルの読み込が成功
	///		  !=0 読み込みに失敗した(std::ifstreamのエラーコード)
	/// 事前条件:既にオープンされたファイルストリームを渡す.
	/// 事後条件:vecPosに全ての座標値が格納され、 m_boxにボックス枠の座標値が設定される.
	///			ファイルストリームが進む
	int LoadBinary(
		std::ifstream& in, ///< [in/out]:既にオープンされたファイルストリーム.
		std::vector<MGPosition>& vecPos ///<[out]:ファイルから読み込んだ座標値の配列.
		);
	
	/// ファイルから読み込んだ全ての頂点の座標値が格納されている.
	/// vecPosから各三角形の法線を計算し、m_vecNormlTriangに追加.
	/// vecPosの各要素の座標値の重複をトレランスを元に判断し
	/// 重複を取り除き、m_vecPosに追加.
	/// vecPosの各要素のインデックスをm_indicesに追加.
	/// 事後条件:m_vecNormlTriang, m_indicesに値が設定される.
	void set_mesh_data(
		const std::vector<MGPosition>& vecPos ///< [in]:ファイルから読み込んだ座標値の配列.
	);

	/// 入力であるpositionがVertexMapにすでに登録されているかチェックする.
	/// 登録されていればそのm_vecPos添え字を返し,
	/// 未登録の場合、新たにm_vecPosに格納し、positionとm_vecPos添え字のmapをVertexMapに
	/// 登録する.
	int IdentifyPosition(
		const MGPosition& position, ///< [in]:頂点の座標.
		triangleMap& VertexMap ///< [in/out]:頂点の座標、頂点のインデックスを保持する.
	);

	///TLDataの読み込みを行う.
	void AddTLData(const mgTLData& tlData, triangleMap& VertexMap);

	//////////

	/// STLファイルの図形のbounding boxを格納する.
	MGBox m_box;

	/// STLファイルの図形を構成する各三角形の、座標の重複がない頂点座標の配列.
	/// ファイルから読み込んだ順で、座標の重複を取り除き、座標値が格納されている.
	std::vector<MGPosition> m_vecPos;

	/// STLファイルの図形を構成する各三角形の法線ベクトルの配列.
	/// ファイルから読み込まれた順で三角形の法線ベクトルが格納されている.
	///m_vecNormlTriang.size()*3=m_indices.size().
	///m_vecNormlTriang[i] is the normal of the triangle  m_indices[i], [i+1], [i+2] for
	///i=0, ..., m_indices.size()/3.
	std::vector<MGUnit_vector> m_vecNormlTriang;

	/// STLファイルの図形を構成する各三角形の各頂点に対応する
	/// 座標の配列のインデックスを格納する配列.
	/// ファイルから読み込んだ順番で各三角形の頂点の並びが格納されている.
	/// 各要素には該当する頂点座標の配列の添え字が格納されている.
	/// 例：i番目の三角形の各頂点の座標の配列の添え字は
	/// m_indices[i*3]、[i*3+1]、[i*3+2]という以上の順序で取得できる.
	std::vector<int> m_indices;
};

/** @} */ // end of MGObjectRelated group

#endif
